# %load run_single_band.py
from triqs.gf import *
from triqs.operators import *
from h5 import *
from triqs_cthyb import Solver
import triqs.utility.mpi as mpi
import numpy as np
import triqs.utility.mpi as mpi
from h5 import HDFArchive
from triqs.operators import *
from triqs_cthyb import *
from triqs.stat.histograms import Histogram
from triqs.gf import *
import numpy as np
from triqs.stat.histograms import *
from triqs.plot.mpl_interface import plt, oplot
from matplotlib.backends.backend_pdf import PdfPages

import os

# Parameters
D, V, U = 1.0, 0.7, 8.0
e_f, beta = -4.0, 4
#e_f, beta = -U/2.0, 10

# Construct the impurity solver with the inverse temperature
# and the structure of the Green's functions
S = Solver(beta = beta, gf_struct = [ ('up',[0]), ('down',[0]) ])

# Initialize the non-interacting Green's function S.G0_iw
for name, g0 in S.G0_iw: g0 << inverse(iOmega_n - e_f - V**2 * Wilson(D))

# Run the solver. The results will be in S.G_tau, S.G_iw and S.G_l
S.solve(h_int = U * n('up',0) * n('down',0),   # Local Hamiltonian 
n_cycles  = 50000000,                           # Number of QMC cycles
length_cycle=50,
move_double=False,
move_shift=False,
wang_landau_cycle  = True, 
wang_landau_lambda  = 1.1000,
n_warmup_cycles = 200000,                      # Warmup cycles
performance_analysis = True,
measure_G2_tau=False,
)

arch = HDFArchive('singleband_hubbard.h5','w')
arch.create_group('CTHYB')
gr = arch['CTHYB']
gr['G_tau'] = S.G_tau
gr['G_tau_accum'] = S.G_tau_accum
gr['asymmetry_G_tau'] = S.asymmetry_G_tau
gr['beta'] = beta
gr['U'] = U
#gr['perturbation_order'] = S.perturbation_order
#gr['perturbation_order_total'] = S.perturbation_order_total
gr['performance_analysis'] = S.performance_analysis


pp = PdfPages('G_asymm_bath.pdf')
sumed1=0;sumed2=0;sumed3=0;sumed4=0
for e_group_name in arch:
    e_group = arch[e_group_name]

    def plot_histos(title,histos,sumed):
        a = plt.gca()
        a.set_xlim((0,beta))
        a.set_xlabel('$\\Delta\\tau$')
        a.set_ylabel(title)

        #sumed=0
        for name, histo in histos.items():
            w = histo.data
            print(name, sum(w))
            sumed = sumed+sum(w)
            dtau = [histo.mesh_point(n) for n in range(len(histo))]
            plt.plot(dtau,w,label=name,linewidth=0.7)
            a.legend(loc='upper center',prop={'size':10})
#        print(sumed)


    histo = e_group['performance_analysis']
    print("Insertion statistics: ")
    plt.subplot(4,1,1)
    proposed = histo['insert_length_proposed_up'] + histo['insert_length_proposed_down']
    accepted = histo['insert_length_accepted_up'] + histo['insert_length_accepted_down']
    plot_histos("Insertion",{"Proposed" : proposed, "Accepted" : accepted},sumed1)
    # Move remove
    print("Removal statistics: ")
    plt.subplot(4,1,2)
    proposed = histo['remove_length_proposed_up'] + histo['remove_length_proposed_down']
    accepted = histo['remove_length_accepted_up'] + histo['remove_length_accepted_down']
    plot_histos("Removal",{"Proposed" : proposed, "Accepted" : accepted},sumed2)
    # worm move insert
    print("Worm Insertion statistics: ")
    plt.subplot(4,1,3)
    proposed = histo['worm_insert_length_proposed_up'] + histo['worm_insert_length_proposed_down']
    accepted = histo['worm_insert_length_accepted_up'] + histo['worm_insert_length_accepted_down']
    plot_histos("Worm Insertion",{"Proposed" : proposed, "Accepted" : accepted},sumed3)
    # worm Move remove
    print("Worm Removal statistics: ")
    plt.subplot(4,1,4)
    proposed = histo['worm_remove_length_proposed_up'] + histo['worm_remove_length_proposed_down'] 
    accepted = histo['worm_remove_length_accepted_up'] + histo['worm_remove_length_accepted_down'] 
    plot_histos("Worm Removal",{"Proposed" : proposed, "Accepted" : accepted},sumed4)
#    # double move insert
#    print("Double Insertion statistics: ")
#    plt.subplot(6,1,5)
#    proposed = histo['double_insert_length_proposed_up'] + histo['double_insert_length_proposed_down']
#    accepted = histo['double_insert_length_accepted_up'] + histo['double_insert_length_accepted_down']
#    plot_histos("Double Insertion",{"Proposed" : proposed, "Accepted" : accepted},sumed1)
#    # double Move remove
#    print("Double Removal statistics: ")
#    plt.subplot(6,1,6)
#    proposed = histo['double_remove_length_proposed_up'] + histo['double_remove_length_proposed_down']
#    accepted = histo['double_remove_length_accepted_up'] + histo['double_remove_length_accepted_down']
#    plot_histos("Double Removal",{"Proposed" : proposed, "Accepted" : accepted},sumed2)

#    print("total sumed: ", sumed1+sumed2+sumed3+sumed4)

    pp.savefig(plt.gcf())
pp.close()
