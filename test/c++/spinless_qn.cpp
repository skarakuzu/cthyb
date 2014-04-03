#include "./ctqmc_matrix.hpp"
#include "./operator.hpp"
#include "./fundamental_operator_set.hpp"
#include <triqs/gfs/local/fourier_matsubara.hpp>
#include <triqs/parameters.hpp>
#include <triqs/gfs/block.hpp>
#include <triqs/gfs/imtime.hpp>
#include <triqs/gfs/imfreq.hpp>

using namespace cthyb_matrix;
using triqs::utility::many_body_operator;
using triqs::utility::c;
using triqs::utility::c_dag;
using triqs::utility::n;
using triqs::utility::parameters;
using namespace triqs::gfs;

int main(int argc, char* argv[]) {

  std::cout << "Welcome to the CTHYB solver\n";

  // Initialize mpi
  boost::mpi::environment env(argc, argv);
  int rank;
  {
    boost::mpi::communicator c;
    rank = c.rank();
  }

  // Parameters
  double beta = 10.0;
  double U = 2.0;
  double mu = 1.0;
  double epsilon = 2.3;
  double t = 0.1;

  // Put in the class
  parameters p;
  p["beta"] = beta;
  p["max_time"] = -1;
  p["random_name"] = "";
  p["random_seed"] = 123 * rank + 567;
  p["max_time"] = -1;
  p["verbosity"] = 3;
  p["length_cycle"] = 50;
  p["n_warmup_cycles"] = 50;
  p["n_cycles"] = 3000;
  p["n_tau_delta"] = 1000;
  p["n_tau_g"] = 1000;

  // define operators
  auto H = U*n("0")*n("1") - mu*(n("0")+n("1")) - t*c_dag("0")*c("1") - t*c_dag("1")*c("0");

  // quantum numbers
  std::vector<many_body_operator<double>> qn;
  qn.push_back(n("0")+n("1"));

  // basis of operators to use
  fundamental_operator_set fops;
  fops.insert("0");
  fops.insert("1");

  // block structure of GF
  std::vector<block_desc_t> block_structure;
  block_structure.push_back({"tot", {{"0"}, {"1"}}});

  // Construct CTQMC solver
  ctqmc_matrix solver(p, H, qn, fops, block_structure);

  // Set hybridization function
  triqs::clef::placeholder<0> om_;
  auto delta_w = gf<imfreq>{{beta, Fermion}, {2,2}};
  auto d00 = slice_target(delta_w, triqs::arrays::range(0,1), triqs::arrays::range(0,1));
  auto d11 = slice_target(delta_w, triqs::arrays::range(1,2), triqs::arrays::range(1,2));
  auto d01 = slice_target(delta_w, triqs::arrays::range(0,1), triqs::arrays::range(1,2));
  auto d10 = slice_target(delta_w, triqs::arrays::range(1,2), triqs::arrays::range(0,1));
  d00(om_) << (om_-epsilon)*(1.0/(om_-epsilon-t))*(1.0/(om_-epsilon+t)) +(om_+epsilon)*(1.0/(om_+epsilon-t))*(1.0/(om_+epsilon+t));
  d11(om_) << (om_-epsilon)*(1.0/(om_-epsilon-t))*(1.0/(om_-epsilon+t)) +(om_+epsilon)*(1.0/(om_+epsilon-t))*(1.0/(om_+epsilon+t));
  d01(om_) << -t*(1.0/(om_-epsilon-t))*(1.0/(om_-epsilon+t)) -t*(1.0/(om_+epsilon-t))*(1.0/(om_+epsilon+t));
  d10(om_) << -t*(1.0/(om_-epsilon-t))*(1.0/(om_-epsilon+t)) -t*(1.0/(om_+epsilon-t))*(1.0/(om_+epsilon+t));
  solver.deltat_view()[0] = inverse_fourier(delta_w);
  
  // Solve!
  solver.solve(p);
  
  // Save the results
  if(rank==0){
    H5::H5File G_file("spinless_qn.output.h5",H5F_ACC_TRUNC);
    h5_write(G_file,"G_tau",solver.gt_view()[0]);
  }

  return 0;
}