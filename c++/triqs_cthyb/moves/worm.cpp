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

#include "./worm.hpp"

namespace triqs_cthyb {

  histogram *move_worm::add_histo(std::string const &name, histo_map_t *histos) {
    if (!histos) return nullptr;
    auto new_histo = histos->insert({name, {.0, config.beta(), 100}});
    return &(new_histo.first->second);
  }

  move_worm::move_worm(int block_index, int block_size, std::string const &block_name, qmc_data &data, wl_data& data_wl, bool yes_worm, mc_tools::random_generator &rng, histo_map_t *histos)
  : data(data),
    data_wl(data_wl),
    config(data.config),
    rng(rng),
    worm_insert(block_index, block_size, block_name, data, data_wl, yes_worm, rng, histos),
    worm_remove(block_index, block_size, block_name, data, data_wl, yes_worm, rng, histos),
    histo_proposed_worm_insert(add_histo("worm_insert_length_proposed_" + block_name, histos)),
    histo_accepted_worm_insert(add_histo("worm_insert_length_accepted_" + block_name, histos)),
    histo_proposed_worm_remove(add_histo("worm_remove_length_proposed_" + block_name, histos)),
    histo_accepted_worm_remove(add_histo("worm_remove_length_accepted_" + block_name, histos))
  {}
  mc_weight_t move_worm::attempt() {
 

//    std::cout<<"Attempting a worm move !!!, is worm container empty?"<<data_wl.no_worm()<<std::endl;
    //if(data_wl.no_worm())
    if(data_wl.chose_insert_worm())
    {
    attempt_was_insert = true;
    //std::cout<<"Attempting a worm insertion !!!"<<std::endl;
    auto res =  worm_insert.attempt();
    
    if (histo_proposed_worm_insert) *histo_proposed_worm_insert <<  worm_insert.get_dtau();
    
//    data_wl.update_current_space();

    //return worm_insert.attempt();
     //return res * data_wl.get_mu_space();
     return res;

    }
    else
    {
    attempt_was_insert = false;
    //std::cout<<"Attempting a worm remove !!!"<<std::endl;
    
    auto res =  worm_remove.attempt();

    if (histo_proposed_worm_remove) *histo_proposed_worm_remove <<  worm_remove.get_dtau();
//    std::cout<<"histo_proposed_worm_remove: "<<worm_remove.get_dtau()<<std::endl;

    //return worm_remove.attempt();
     //return res / data_wl.get_mu_space();
     return res;
    }

  }

  mc_weight_t move_worm::accept() {

    if(attempt_was_insert)
    {

    auto res =  worm_insert.accept();
    
    if (histo_accepted_worm_insert) *histo_accepted_worm_insert << worm_insert.get_dtau();
    
//    data_wl.update_mu_space();

    return res;
    //return worm_insert.accept();
    }
    else
    {
    
    auto res =  worm_remove.accept();
    
    if (histo_accepted_worm_remove ) *histo_accepted_worm_remove << worm_remove.get_dtau();
//    std::cout<<"histo_accepted_worm_remove: "<<worm_remove.get_dtau()<<std::endl;

//    data_wl.update_current_space();
//    data_wl.update_mu_space();

    return res;
    //return worm_remove.accept();

    }
  }

  void move_worm::reject() {
  
    if(attempt_was_insert)
    {

    worm_insert.reject();

//    data_wl.update_current_space();
//    data_wl.update_mu_space();
    
    }
    else
    {
    
    worm_remove.reject();
//    data_wl.update_mu_space();

    }
  }
}
