/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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

#include "./replace.hpp"

namespace triqs_cthyb {

  histogram *move_worm_replace_operator::add_histo(std::string const &name, histo_map_t *histos) {
    if (!histos) return nullptr;
    auto new_histo = histos->insert({name, {.0, config.beta(), 100}});
    return &(new_histo.first->second);
  }

  move_worm_replace_operator::move_worm_replace_operator(qmc_data &data, wl_data& data_wl, mc_tools::random_generator &rng, histo_map_t *histos)
     : data(data),
       data_wl(data_wl),
       config(data.config),
       rng(rng),
       histo_proposed(add_histo("replace_worm_length_proposed", histos)),
       histo_accepted(add_histo("replace_worm_length_accepted", histos)),
       block_index(0) {}

  mc_weight_t move_worm_replace_operator::attempt() {

#ifdef EXT_DEBUG
    std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cerr << "In config " << config.get_id() << std::endl;
    std::cerr << "* Attempt for move_shift_operator ";
#endif

    // --- Choose an operator in configuration to shift at random
    // By choosing an *operator* in config directly, not bias based on det size introduced
    auto worm_size = data_wl.size_worm();
    if (worm_size == 0) {
#ifdef EXT_DEBUG
      std::cerr << "(empty configuration)" << std::endl;
      block_index = -1;
#endif
      return 0;
    }
    const int op_pos_in_worm = rng(worm_size);

    auto is_dagger = rng(2) == 1 ? true : false;

    // --- Find operator (and its characteristics) from the configuration
    tau_worm        = is_dagger ? data_wl.get_time_dag(op_pos_in_worm) : data_wl.get_time(op_pos_in_worm);
    op_worm         = is_dagger ? data_wl.get_op_dag(op_pos_in_worm) : data_wl.get_op(op_pos_in_worm);
    block_index    = op_worm.block_index;

#ifdef EXT_DEBUG
    std::cerr << "(block " << block_index << ")" << std::endl;
#endif

    // Properties corresponding to det
    auto &det     = data.dets[block_index];
    auto det_size = det.size();
    if (det_size == 0) return 0; // nothing to replace

    // Construct new operator
    // Choose a new inner index (this is done here for compatibility)
    auto inner_new = rng(data.n_inner[block_index]);
    op_det         = op_desc{block_index, inner_new, is_dagger, data.linindex[std::make_pair(block_index, inner_new)]};


    const int op_pos_in_det = rng(det_size);


#ifdef EXT_DEBUG
    std::cerr << "* Proposing to replace worm:" << std::endl;
    std::cerr << op_worm << " tau = " << tau_worm << std::endl;
    std::cerr << " to " << std::endl;
    std::cerr << op_det << " tau = " << tau_det << std::endl;
#endif
    std::cout << "* Proposing to replace worm:" << std::endl;
    std::cout << op_worm << " tau = " << tau_worm << std::endl;
    std::cout << " to " << std::endl;
    std::cout << op_det << " tau = " << tau_det << std::endl;

    std::cout<<"printing full configuration....."<<std::endl;
    std::cout<<config<<std::endl;
    // --- Compute the det ratio
/*
    // Do we need to roll the determinant?
    roll_direction = det_type::None;

    // Check if we went through \tau = 0 or \tau = \beta
    if (det_size > 1) {
      if ((tau_old > tL) && (tau_new < tL)) roll_direction = (is_dagger ? det_type::Up : det_type::Left);
      if ((tau_old < tL) && (tau_new > tL)) roll_direction = (is_dagger ? det_type::Down : det_type::Right);
    }
*/
    // Replace old row/column with new operator time/inner_index. Returns the ratio of dets (Cf det_manip doc).
    auto det_ratio = (is_dagger ? det.try_change_row(op_pos_in_det, {tau_worm, op_worm.inner_index}) :
                                  det.try_change_col(op_pos_in_det, {tau_worm, op_worm.inner_index}));

    // for quick abandon
    double random_number = rng.preview();
    if (random_number == 0.0) return 0;

    /*
    double p_yee = std::abs(det_ratio / data.atomic_weight);
    // --- Compute the atomic_weight ratio
    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
    if (new_atomic_weight == 0.0) {
#ifdef EXT_DEBUG
      std::cerr << "atomic_weight == 0" << std::endl;
#endif
      return 0;
    }
    auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "atomic_weight_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " "
                          << new_atomic_weight / data.atomic_weight << " in config " << config.get_id();
  */

    // --- Compute the weight
    mc_weight_t p = det_ratio;

#ifdef EXT_DEBUG
    std::cerr << "Det ratio: " << det_ratio << '\t';
#endif

    return p;

  }

  mc_weight_t move_worm_replace_operator::accept() {

    /*
    // Update the tree
    data.imp_trace.confirm_shift();

    // Update the configuration
    config.erase(tau_old);
    config.insert(tau_new, op_new);
    config.finalize();

    // Update the determinant
    data.dets[block_index].complete_operation();
    data.update_sign();

    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;

    if (histo_accepted) *histo_accepted << dtau;

    auto result = data.current_sign / data.old_sign * data.dets[block_index].roll_matrix(roll_direction);

#ifdef EXT_DEBUG
    std::cerr << "* Move move_shift_operator accepted" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    check_det_sequence(data.dets[block_index], config.get_id());
#endif
*/
    return mc_weight_t(1.);
  }

  void move_worm_replace_operator::reject() {
/*
    config.finalize();
    data.imp_trace.cancel_shift();
    data.dets[block_index].reject_last_try();
*/
#ifdef EXT_DEBUG
    std::cerr << "* Move move_shift_operator rejected" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    if (block_index != -1) check_det_sequence(data.dets[block_index], config.get_id());
#endif
  }
}
