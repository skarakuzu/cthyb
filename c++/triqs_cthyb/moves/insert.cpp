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

#include "./insert.hpp"

namespace triqs_cthyb {

  histogram *move_insert_c_cdag::add_histo(std::string const &name, histo_map_t *histos) {
    if (!histos) return nullptr;
    auto new_histo = histos->insert({name, {.0, config.beta(), 100}});
    return &(new_histo.first->second);
  }

  move_insert_c_cdag::move_insert_c_cdag(int block_index, int block_size, std::string const &block_name, qmc_data &data, wl_data& data_wl, bool yes_worm, mc_tools::random_generator &rng, histo_map_t *histos)
     : data(data),
       data_wl(data_wl),
       yes_worm(yes_worm),
       config(data.config),
       rng(rng),
       block_index(block_index),
       block_size(block_size),
       histo_proposed(add_histo("insert_length_proposed_" + block_name, histos)),
       histo_accepted(add_histo("insert_length_accepted_" + block_name, histos)) {}

  mc_weight_t move_insert_c_cdag::attempt() {


   /*
    std::cout<<std::endl;
    if(!yes_worm) std::cout<<"************NEW INSERT MOVE ATTEMPT********"<<std::endl;
    else std::cout<<"************NEW WORM INSERT MOVE ATTEMPT********"<<std::endl;
    std::cout<<std::endl;
    */

#ifdef EXT_DEBUG
    std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cerr << "In config " << config.get_id() << std::endl;
    std::cerr << "* Attempt for move_insert_c_cdag (block " << block_index << ")" << std::endl;
#endif

    // Pick up the value of alpha and choose the operators
    auto rs1 = rng(block_size), rs2 = rng(block_size);
    op1 = op_desc{block_index, rs1, true, data.linindex[std::make_pair(block_index, rs1)]};
    op2 = op_desc{block_index, rs2, false, data.linindex[std::make_pair(block_index, rs2)]};

    // Choice of times for insertion. Find the time as double and them put them on the grid.
    tau1 = data.tau_seg.get_random_pt(rng);
    tau2 = data.tau_seg.get_random_pt(rng);

#ifdef EXT_DEBUG
    std::cerr << "* Proposing to insert:" << std::endl;
    std::cerr << op1 << " at " << tau1 << std::endl;
    std::cerr << op2 << " at " << tau2 << std::endl;
#endif
    /*
    std::cout << "* Proposing to insert:" << std::endl;
    std::cout << op1 << " at " << tau1 << std::endl;
    std::cout << op2 << " at " << tau2 << std::endl;
    */

    // record the length of the proposed insertion
    dtau = double(tau2 - tau1);
    if (histo_proposed && !yes_worm) *histo_proposed << dtau;

    // Insert the operators op1 and op2 at time tau1, tau2
    // 1- In the very exceptional case where the insert has failed because an operator is already sitting here
    // (cf std::map doc for insert return), we reject the move.
    // 2- If ok, we store the iterator to the inserted operators for later removal in reject if necessary
    try {
      data.imp_trace.try_insert(tau1, op1);
      data.imp_trace.try_insert(tau2, op2);
    } catch (rbt_insert_error const &) {
      std::cerr << "Insert error : recovering ... " << std::endl;
      data.imp_trace.cancel_insert();
      return 0;
    }

    
    
    int det_size=0;
    mc_weight_t t_ratio, det_ratio;
    double p_yee, random_number;
    
    /*
    std::cout<<"printing full configuration....."<<yes_worm<<std::endl;
    std::cout<<config<<std::endl;
    */


    if(!yes_worm)
    {
    // Computation of det ratio
    auto &det    = data.dets[block_index];
    int det_size = det.size();

    // Find the position for insertion in the determinant
    // NB : the determinant stores the C in decreasing time order.
    int num_c_dag, num_c;
    for (num_c_dag = 0; num_c_dag < det_size; ++num_c_dag) {
      if (det.get_x(num_c_dag).first < tau1) break;
    }
    for (num_c = 0; num_c < det_size; ++num_c) {
      if (det.get_y(num_c).first < tau2) break;
    }
/*   
   if(num_c!=0 && num_c_dag!=0)
   { 
    std::cerr << "* testing det functions inside insert:" << std::endl;
    std::cerr << "num_c_dag " << " is " << num_c_dag<< std::endl;
    std::cerr << "num_c " << " is " << num_c << std::endl;
    //std::cerr << "num_c_dag vs det.get_x(num_c_dag).first" << " is " << num_c_dag<<" "<< det.get_x(num_c_dag).first << std::endl;
    //std::cerr << "num_c vs det.get_y(num_c).first" << " is " << num_c<<" "<< det.get_y(num_c).first << std::endl;
	}
*/
    // Insert in the det. Returns the ratio of dets (Cf det_manip doc).
    det_ratio = det.try_insert(num_c_dag, num_c, {tau1, op1.inner_index}, {tau2, op2.inner_index});
    
    // proposition probability
    t_ratio = std::pow(block_size * config.beta() / double(det.size() + 1), 2);

    // For quick abandon
    random_number = rng.preview();
    if (random_number == 0.0) return 0;
    p_yee = std::abs(t_ratio * det_ratio / data.atomic_weight);

     
    }
    else
    {

//      std::cout<<"Inside worm insertion in insert.cpp"<<std::endl;
    det_ratio= mc_weight_t(1.0);
    
    data_wl.insert_worm_dag(tau1, op1);
    data_wl.insert_worm(tau2, op2);
    
    data_wl.update_current_space();

    // proposition probability
    t_ratio = std::pow(block_size * config.beta() , 2);

    // For quick abandon
    random_number = rng.preview();
    if (random_number == 0.0) return 0;
    p_yee = std::abs(t_ratio * det_ratio / data.atomic_weight);

    }

    p_yee = -1;

    // computation of the new trace after insertion
    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
    if (new_atomic_weight == 0.0) {
#ifdef EXT_DEBUG
      std::cerr << "atomic_weight == 0" << std::endl;
#endif
 //     std::cout << "atomic_weight == 0 with pyee " << p_yee<<" "<<random_number<< std::endl;
      return 0;
    }
    auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "(insert) trace_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " "
                          << new_atomic_weight / data.atomic_weight << " in config " << config.get_id();

    
     mc_weight_t p = atomic_weight_ratio * det_ratio;

#ifdef EXT_DEBUG
    std::cerr << "Atomic ratio: " << atomic_weight_ratio << '\t';
    std::cerr << "Det ratio: " << det_ratio << '\t';
    std::cerr << "Prefactor: " << t_ratio << '\t';
    std::cerr << "Weight: " << p * t_ratio << std::endl;
    std::cerr << "p_yee * newtrace: " << p_yee * new_atomic_weight << std::endl;
#endif

    if (!isfinite(p * t_ratio)) {
      std::cerr << "Insert move info:\n";
      std::cerr << "Atomic ratio: " << atomic_weight_ratio << '\t';
      if(!yes_worm) std::cerr << "Det ratio: " << det_ratio << '\t';
      std::cerr << "Prefactor: " << t_ratio << '\t';
      std::cerr << "Weight: " << p * t_ratio << std::endl;
      std::cerr << "p_yee * newtrace: " << p_yee * new_atomic_weight << std::endl;
      
      TRIQS_RUNTIME_ERROR << "(insert) p * t_ratio not finite p : " << p << " t_ratio : " << t_ratio << " in config " << config.get_id();
    }

    //just to test indertion; delete later
     //if(yes_worm) return mc_weight_t(1.0);
     //if(yes_worm) std::cout<<"worm insertion ratio***** : "<<p*t_ratio<<std::endl;
     //else std::cout<<" insertion ratio***** : "<<p*t_ratio<<std::endl;
   


     if(yes_worm) return p * t_ratio * data_wl.get_mu_space();
     else return p * t_ratio ;
     //return p * t_ratio * data_wl.get_mu_space();
  }

  mc_weight_t move_insert_c_cdag::accept() {

    // insert in the tree
    data.imp_trace.confirm_insert();

    // insert in the configuration
    config.insert(tau1, op1);
    config.insert(tau2, op2);
    config.finalize();
/*
    if(yes_worm)
    {
      std::cout<<"HERE inserting accepted worm"<<std::endl;
    data_wl.insert_worm_dag(tau1, op1);
    data_wl.insert_worm(tau2, op2);
    }
*/
    data_wl.update_current_space();
    data_wl.update_mu_space();
    data_wl.update_visit_space();
    
    // insert in the determinant
    if(!yes_worm) data.dets[block_index].complete_operation();
    //data.update_sign();
    if(data_wl.no_worm()) data.update_sign();
    else data.update_my_sign(data_wl);
    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;
    if (histo_accepted && !yes_worm) *histo_accepted << dtau;

#ifdef EXT_DEBUG
    std::cerr << "* Move move_insert_c_cdag accepted" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    check_det_sequence(data.dets[block_index], config.get_id());
#endif
//    check_det_sequence(data.dets[block_index], config.get_id());

    return data.current_sign / data.old_sign;
  }

  void move_insert_c_cdag::reject() {

    config.finalize();
    data.imp_trace.cancel_insert();
    if(!yes_worm) data.dets[block_index].reject_last_try();
    
    if(yes_worm && !data_wl.no_worm())
    {
      //std::cout<<"HERE deleting rejected worm"<<std::endl;
    data_wl.erase_worm_dag(0);
    data_wl.erase_worm(0);

    }
    data_wl.update_current_space();
    data_wl.update_mu_space();
    data_wl.update_visit_space();

#ifdef EXT_DEBUG
    std::cerr << "* Move move_insert_c_cdag rejected" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    check_det_sequence(data.dets[block_index], config.get_id());
#endif
//    check_det_sequence(data.dets[block_index], config.get_id());
  }
}
