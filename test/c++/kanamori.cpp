#include "./ctqmc_krylov.hpp"
#include "./operator.hpp"
#include "./fundamental_operator_set.hpp"
#include <triqs/gfs/local/fourier_matsubara.hpp>
#include <triqs/parameters.hpp>
#include <triqs/gfs/block.hpp>
#include <triqs/gfs/imtime.hpp>
#include <triqs/gfs/imfreq.hpp>

using namespace triqs::app::impurity_solvers::ctqmc_krylov;
using triqs::utility::many_body_operator;
using triqs::utility::c;
using triqs::utility::c_dag;
using triqs::utility::n;
using triqs::utility::parameters;
using triqs::gfs::gf;
using triqs::gfs::imfreq;
using triqs::gfs::imfreq;
using triqs::gfs::make_gf;
using triqs::gfs::Fermion;

int main(int argc, char* argv[]) {

  std::cout << "Welcome to the Krylov solver\n";

  // Initialize mpi
  boost::mpi::environment env(argc, argv);
  int rank;
  {
    boost::mpi::communicator c;
    rank = c.rank();
  }

  // Parameters
  double beta = 10.0;
  int num_orbitals = 2;
  double mu = 1.0;
  double U = 2.0;
  double J = 0.2;
  double V = 1.0;
  double epsilon = 2.3;

  // Put in the class
  parameters p;
  p["beta"] = beta;
  p["Max_Time"] = -1;
  p["Random_Generator_Name"] = "";
  p["Random_Seed"] = 123 * rank + 567;
  p["Verbosity"] = 3;
  p["Length_Cycle"] = 50;
  p["N_Warmup_Cycles"] = 50;
  p["N_Cycles"] = 500;
  p["n_tau_delta"] = 1000;
  p["n_tau_g"] = 1000;
  p["krylov_bs_use_cutoff"] = true;
  p["krylov_bs_prob_cutoff"] = .0;
  
  // basis of operators to use  
  fundamental_operator_set<int, const char *> fops;
  for(int o = 0; o < num_orbitals; ++o){
      fops.add_operator(o,"up");
      fops.add_operator(o,"down");
  }

  // block structure of GF
  std::vector<block_desc_t<int,const char *>> block_structure;
  for(int o = 0; o < num_orbitals; ++o){
    std::stringstream bup; bup << "up-" << o;
    block_structure.push_back({bup.str(),{std::make_tuple(o,"up")}});
  }
  for(int o = 0; o < num_orbitals; ++o){
    std::stringstream bdown; bdown << "down-" << o;
    block_structure.push_back({bdown.str(),{std::make_tuple(o,"down")}});
  }
    
  // Hamiltonian
  many_body_operator<double,int,const char*> H;
  for(int o = 0; o < num_orbitals; ++o){
      H += -mu*(n(o,"up") + n(o,"down"));
  }
  for(int o = 0; o < num_orbitals; ++o){
      H += U *n(o,"up")*n(o,"down");
  }
  for(int o1 = 0; o1 < num_orbitals; ++o1)
  for(int o2 = 0; o2 < num_orbitals; ++o2){
      if(o1==o2) continue;
      H += (U-2*J)*n(o1,"up")*n(o2,"down");
  }
  for(int o1 = 0; o1 < num_orbitals; ++o1)
  for(int o2 = 0; o2 < num_orbitals; ++o2){
      if(o2>=o1) continue;
      H += (U-3*J)*n(o1,"up")*n(o2,"up");
      H += (U-3*J)*n(o1,"down")*n(o2,"down");
  }
  
  for(int o1 = 0; o1 < num_orbitals; ++o1)
  for(int o2 = 0; o2 < num_orbitals; ++o2){
      if(o1==o2) continue;
      H += -J*c_dag(o1,"up")*c_dag(o1,"down")*c(o2,"up")*c(o2,"down");
      H += -J*c_dag(o1,"up")*c_dag(o2,"down")*c(o2,"up")*c(o1,"down");
  }

  // quantum numbers
  std::vector<many_body_operator<double,int,const char*>> qn;

  // Construct CTQMC solver
  ctqmc_krylov solver(p, H, qn, fops, block_structure);

  // Set hybridization function
  triqs::clef::placeholder<0> om_;
  auto delta_w = make_gf<imfreq>(beta, Fermion, make_shape(1,1));
  delta_w(om_) << V / (om_ - epsilon) + V / (om_ + epsilon);  
  for (int o = 0; o < 2*num_orbitals; ++o){
    solver.deltat_view()[o] = triqs::gfs::lazy_inverse_fourier(delta_w);
  }

  // Solve!
  solver.solve(p);
  
  // Save the results
  if(rank==0){
    H5::H5File G_file("kanamori.output.h5",H5F_ACC_TRUNC);
    for(int o = 0; o < num_orbitals; ++o) {
      std::stringstream bup; bup << "G_up-" << o;
      h5_write(G_file, bup.str(), solver.gt_view()[o]);
      std::stringstream bdown; bdown << "G_down-" << o;
      h5_write(G_file, bdown.str(), solver.gt_view()[num_orbitals+o]);
    }
  }

  return 0;

}
