/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include "ctqmc_matrix.hpp"
#include <triqs/utility/callbacks.hpp>

#include "move_insert.hpp"
#include "move_remove.hpp"
#include "measure_g.hpp"
#include "measure_perturbation_hist.hpp"

namespace cthyb_matrix {

ctqmc_matrix::ctqmc_matrix(parameters p_in, real_operator_t const& h_loc, std::vector<real_operator_t> const& quantum_numbers,
                           fundamental_operator_set const& fops, std::vector<block_desc_t> const& block_structure)
   : gf_block_structure(fops, block_structure) {
 p_in.update(constructor_defaults());//, utility::parameters::reject_key_without_default);
 auto const& params = p_in;

 if (params["use_quantum_numbers"]) 
  sosp = {h_loc, quantum_numbers, fops};
 else 
  sosp = {h_loc, fops};

 std::vector<std::string> block_names;
 std::vector<gf<imtime>> deltat_blocks;
 std::vector<gf<imtime>> gt_blocks;

 for (auto const& block : block_structure) {
  block_names.push_back(block.name);
  int n = block.indices.size();
  deltat_blocks.push_back(gf<imtime>{{params["beta"], Fermion, params["n_tau_delta"], half_bins}, {n, n}});
  gt_blocks.push_back(gf<imtime>{{params["beta"], Fermion, params["n_tau_g"], half_bins}, {n, n}});
 }

 deltat = make_block_gf(block_names, deltat_blocks);
 gt = make_block_gf(block_names, gt_blocks);
}

//-----------------------------------

void ctqmc_matrix::solve(utility::parameters p_in) {

 p_in.update(solve_defaults());//, utility::parameters::reject_key_without_default);
 auto const& params = p_in;

 qmc_data data(params, sosp, gf_block_structure, deltat);
 mc_tools::mc_generic<mc_sign_type> qmc(params);

 // Moves
 auto& delta_names = deltat.domain().names();
 for (size_t block = 0; block < deltat.domain().size(); ++block) {
  int block_size = deltat[block].data().shape()[1];
  qmc.add_move(move_insert_c_cdag(block, block_size, data, qmc.rng(), false), "Insert Delta_" + delta_names[block]);
  qmc.add_move(move_remove_c_cdag(block, block_size, data, qmc.rng()), "Remove Delta_" + delta_names[block]);
 }

 // Measurements
 if (params["measure_gt"]) {
  auto& gt_names = gt.domain().names();
  for (size_t block = 0; block < gt.domain().size(); ++block) {
   qmc.add_measure(measure_g(block, gt[block], data), "G measure (" + gt_names[block] + ")");
  }
 }
 if (params["measure_pert_order"]) {
  auto& gt_names = gt.domain().names();
  for (size_t block = 0; block < gt.domain().size(); ++block) {
   qmc.add_measure(measure_perturbation_hist(block, data, "histo_pert_order_" + gt_names[block] + ".dat"), "Perturbation order (" + gt_names[block] + ")");
  }
 }

 // run!! The empty configuration has sign = 1
 qmc.start(1.0, triqs::utility::clock_callback(params["max_time"]));
 qmc.collect_results(c);
}

//----------------------------------------------------------------------------------------------
parameter_defaults ctqmc_matrix::constructor_defaults() const {

 parameter_defaults pdef;

 pdef.required("beta", double(), "Inverse temperature")
     .optional("n_tau_delta", int(10001), "Number of time slices for Delta(tau)")
     .optional("n_tau_g", int(10001), "Number of time slices for G(tau)")
     .optional("n_w", int(1025), "Number of Matsubara frequencies")
     .optional("use_quantum_numbers", bool(false), " Use the quantum numbers");
 return pdef;
}

//----------------------------------------------------------------------------------------------

parameter_defaults ctqmc_matrix::solve_defaults() const {

 parameter_defaults pdef;

 pdef.required("n_cycles", int(), "Number of QMC cycles")
     .optional("length_cycle", int(50), "Length of a single QMC cycle")
     .optional("n_warmup_cycles", int(5000), "Number of cycles for thermalization")
     .optional("random_seed", int(34788 + 928374 * c.rank()), "Seed for random number generator")
     .optional("random_name", std::string(""), "Name of random number generator")
     .optional("max_time", int(-1), "Maximum runtime in seconds, use -1 to set infinite")
     .optional("verbosity", (c.rank()==0 ? int(3) : int(0)), "Verbosity level")
     .optional("measure_gt", bool(true), "Whether to measure G(tau)")
     .optional("measure_pert_order", bool(false), "Whether to measure perturbation order")
     .optional("make_histograms", bool(false), "Make the analysis histograms of the trace computation")
     .optional("trace_estimator", std::string("FullTrace"), "How to compute the trace");
 return pdef;
}

void ctqmc_matrix::help() const {
 // TODO
}
}
