import h5py
import numpy as np
#from matplotlib import pyplot as plt
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

U = 8.0
mu = -4.0
beta = 4.0
d = 0.7

#filename = 'singleband_hubbard.h5_nowl'
filename = 'singleband_hubbard.h5_d'+str(d)+'_U'+str(U)+'_bt'+str(beta)+'_mu'+str(mu)+'_nowl'
#filename1 = 'singleband_hubbard.h5_yeswl'
#filename1 = 'singleband_hubbard.h5_d'+str(d)+'_U'+str(U)+'_bt'+str(beta)+'_mu'+str(mu)+'_yeswl_3space'
filename1 = 'singleband_hubbard.h5_d'+str(d)+'_U'+str(U)+'_bt'+str(beta)+'_mu'+str(mu)+'_yeswl_3space'
filename2 = 'singleband_hubbard.h5_d'+str(d)+'_U'+str(U)+'_bt'+str(beta)+'_mu'+str(mu)+'_yeswl_2space'

data = h5py.File(filename,'r')
data1 = h5py.File(filename1,'r')
data2 = h5py.File(filename2,'r')

print('U:',U)
print('mu:',mu)
print('beta:',beta)
print('d:',d)

gf_up_nowl = data['CTHYB/G_tau']['up']['data'][:]
gf_dwn_nowl = data['CTHYB/G_tau']['down']['data'][:]
gf_up_yeswl = data1['CTHYB/G_tau']['up']['data'][:]
gf_dwn_yeswl = data1['CTHYB/G_tau']['down']['data'][:]
gf_up_yeswl_2 = data2['CTHYB/G_tau']['up']['data'][:]
gf_dwn_yeswl_2 = data2['CTHYB/G_tau']['down']['data'][:]



#tau2 = np.linspace(0.0, beta, num=int(Z.shape[0]/2));
tau = np.linspace(0.0, beta, num=10001);
print('tau shape: ',tau.shape)
#print('tau2 shape: ',tau2.shape)

plt.plot(tau,gf_up_nowl[:,0,0,0],'-o', label='CTHYB-nowl')
#plt.plot(tau,gf_dwn_nowl[:,0,0,0],'-o', label='CTHYB-nowl')
plt.plot(tau,gf_up_yeswl[:,0,0,0],'-o', label='CTHYB-yeswl-3space') # Z, Gup, Gdwn
plt.plot(tau,gf_up_yeswl_2[:,0,0,0],'-o', label='CTHYB-yeswl-2space') # Z, G(up+dwn)
#plt.plot(tau,-gf_dwn_yeswl[:,0,0,0],'-o', label='CTHYB-yeswl')


plt.ylim(-0.8, 0.0)
plt.legend(loc="center")

plt.show()
#plt.plot(gf_down[:,0,0,0],'*')
plt.savefig('myplot_gf.png')

#plt.show()


