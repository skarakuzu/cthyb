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

#include "./remove.hpp"

namespace triqs_cthyb {

  histogram * move_remove_c_cdag::add_histo(std::string const &name, histo_map_t *histos) {
    if (!histos) return nullptr;
    auto new_histo = histos->insert({name, {.0, config.beta(), 100}});
    return &(new_histo.first->second);
  }
  
  //move_remove_c_cdag::move_remove_c_cdag(int block_index, int block_size, std::string const &block_name, qmc_data &data, mc_tools::random_generator &rng, histo_map_t *histos)
//     : data(data),

  move_remove_c_cdag::move_remove_c_cdag(int block_index, int block_size, std::string const &block_name, qmc_data &data, wl_data& data_wl, bool yes_worm, mc_tools::random_generator &rng, histo_map_t *histos)
     : data(data),
       data_wl(data_wl),
       yes_worm(yes_worm),
       config(data.config),
       rng(rng),
       block_index(block_index),
       block_size(block_size),
       histo_proposed(add_histo("remove_length_proposed_" + block_name, histos)),
       histo_accepted(add_histo("remove_length_accepted_" + block_name, histos)) {}

  mc_weight_t move_remove_c_cdag::attempt() {

#ifdef EXT_DEBUG
    std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cerr << "In config " << config.get_id() << std::endl;
    std::cerr << "* Attempt for move_remove_c_cdag (block " << block_index << ")" << std::endl;
#endif
	
    int det_size=0;
    mc_weight_t t_ratio, det_ratio;
    double p_yee, random_number;

    if(!yes_worm)
    {

    auto &det = data.dets[block_index];

    // Pick up a couple of C, Cdagger to remove at random
    // Remove the operators from the traces
    det_size = det.size();
    if (det_size == 0) return 0; // nothing to remove
    int num_c_dag = rng(det_size), num_c = rng(det_size);

#ifdef EXT_DEBUG
    std::cerr << "* Proposing to remove: ";
    std::cerr << num_c_dag << "-th Cdag(" << block_index << ",...), ";
    std::cerr << num_c << "-th C(" << block_index << ",...)" << std::endl;
#endif

    // now mark 2 nodes for deletion
    //tau1 = data.imp_trace.try_delete(num_c, block_index, false);
    //tau2 = data.imp_trace.try_delete(num_c_dag, block_index, true);

    tau1 = det.get_y(num_c).first;
    tau2 = det.get_x(num_c_dag).first;
    int n1 = data.imp_trace.try_delete_worm(tau2, block_index, true);
    int n2 = data.imp_trace.try_delete_worm(tau1, block_index, false);
    
    std::cout << "* Proposing to remove: ";
    std::cout << n1 << "-th Cdag(" << block_index << ",...), at time "<<tau2;
    std::cout << n2 << "-th C(" << block_index << ",...), at time " <<tau1<< std::endl;
 
    std::cout<<"printing full configuration....."<<std::endl;
    std::cout<<config<<std::endl;
    
    /*
    tau1_try = det.get_x(num_c_dag).first;
    tau2_try = det.get_y(num_c).first;
    int n1 = data.imp_trace.try_delete_mine(tau2_try, block_index, false);
    int n2 = data.imp_trace.try_delete_mine(tau1_try, block_index, true);

   if(num_c!=0 && num_c_dag!=0)
   { 
    std::cerr << "* testing det functions inside remove:" << std::endl;
    //std::cerr << "num_c_dag " << " is " << num_c_dag<< std::endl;
    //std::cerr << "num_c " << " is " << num_c << std::endl;
    std::cerr << "n2 vs num_c_dag vs det.get_x(num_c_dag).first" << " is " << n2<< " " << num_c_dag<<" "<< det.get_x(num_c_dag).first << std::endl;
    std::cerr << "n1 vs num_c vs det.get_y(num_c).first" << " is " << n1 << " " <<num_c<<" "<< det.get_y(num_c).first << std::endl;
	}
*/


    // record the length of the proposed removal
    dtau = double(tau2 - tau1);
    if (histo_proposed) *histo_proposed << dtau;

    det_ratio = det.try_remove(num_c_dag, num_c);

    // proposition probability
    t_ratio = std::pow(block_size * config.beta() / double(det_size), 2); // Size of the det before the try_delete!

    // For quick abandon
    random_number = rng.preview();
    if (random_number == 0.0) return 0;
    p_yee = std::abs(det_ratio / t_ratio / data.atomic_weight);

    }
    else
    {
   /* my changes */ 
    tau1 = data_wl.get_time_dag(0);
    tau2 = data_wl.get_time(0);
    int block_index_worm_dag = data_wl.get_block_index_dag(0);
    int block_index_worm = data_wl.get_block_index(0);
    
    int num_c_dag = data.imp_trace.try_delete_worm(tau1, block_index_worm_dag, true);
    int num_c = data.imp_trace.try_delete_worm(tau2, block_index_worm, false);
    
    std::cout << "* Proposing to remove worm: ";
    std::cout << num_c_dag << "-th Cdag(" << block_index_worm_dag << ",...), at time "<<tau1;
    std::cout << num_c << "-th C(" << block_index_worm << ",...), at time " <<tau2<< std::endl;

    std::cout<<"printing full configuration....."<<std::endl;
    std::cout<<config<<std::endl;


    
    det_ratio= mc_weight_t(1.0);
    
    t_ratio = std::pow(block_size * config.beta() , 2); // Size of the det before the try_delete!

    // For quick abandon
    random_number = rng.preview();
    if (random_number == 0.0) return 0;
    p_yee = std::abs(det_ratio / t_ratio / data.atomic_weight);
    }


    // recompute the atomic_weight
    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
    if (new_atomic_weight == 0.0) {
#ifdef EXT_DEBUG
      std::cerr << "atomic_weight == 0" << std::endl;
#endif
      return 0;
    }
    auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "(remove) atomic_weight_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " "
                          << new_atomic_weight / data.atomic_weight << " in config " << config.get_id();

    mc_weight_t p = atomic_weight_ratio * det_ratio;

#ifdef EXT_DEBUG
    std::cerr << "Trace ratio: " << atomic_weight_ratio << '\t';
    std::cerr << "Det ratio: " << det_ratio << '\t';
    std::cerr << "Prefactor: " << t_ratio << '\t';
    std::cerr << "Weight: " << p / t_ratio << std::endl;
#endif

    if (!isfinite(p)) {
      std::cerr << "Remove move info\n";
      std::cerr << "Trace ratio: " << atomic_weight_ratio << '\t';
      std::cerr << "Det ratio: " << det_ratio << '\t';
      std::cerr << "Prefactor: " << t_ratio << '\t';
      std::cerr << "Weight: " << p / t_ratio << std::endl;
      TRIQS_RUNTIME_ERROR << "(remove) p not finite :" << p << " in config " << config.get_id();
    }
    
    if (!isfinite(p / t_ratio)){
      TRIQS_RUNTIME_ERROR << "(remove) p / t_ratio not finite p : " << p << " t_ratio :  " << t_ratio << " in config " << config.get_id();
    }
    
    //just to test indertion; delete later
     //if(yes_worm) return mc_weight_t(1.0);
    
     if(yes_worm) std::cout<<"worm removal ratio***** : "<<p/t_ratio<<std::endl;
     else std::cout<<"removal ratio***** : "<<p/t_ratio<<std::endl;
     
     return p / t_ratio;

  }

  mc_weight_t move_remove_c_cdag::accept() {

    // remove from the tree
    data.imp_trace.confirm_delete();

    // remove from the configuration
    config.erase(tau1);
    config.erase(tau2);
    config.finalize();

    if(yes_worm)
    {
      std::cout<<"HERE erasing accepted worm"<<std::endl;
      data_wl.erase_worm(0);
      data_wl.erase_worm_dag(0);
      std::cout<<"worm sizes after removal: "<<data_wl.size_worm()<<" "<<data_wl.size_worm_dag()<<std::endl;
    }

    // remove from the determinants
    if(!yes_worm) data.dets[block_index].complete_operation();
    data.update_sign();
    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;
    if (histo_accepted) *histo_accepted << dtau;

#ifdef EXT_DEBUG
    std::cerr << "* Move move_remove_c_cdag accepted" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    if(!yes_worm) check_det_sequence(data.dets[block_index], config.get_id());
#endif

    return data.current_sign / data.old_sign;
  }

  void move_remove_c_cdag::reject() {

    config.finalize();
    data.imp_trace.cancel_delete();
    if(!yes_worm) data.dets[block_index].reject_last_try();

#ifdef EXT_DEBUG
    std::cerr << "* Move move_remove_c_cdag rejected" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    if(!yes_worm) check_det_sequence(data.dets[block_index], config.get_id());
#endif
  }
}
