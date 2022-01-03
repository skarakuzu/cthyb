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

/*
  histogram *move_worm::add_histo(std::string const &name, histo_map_t *histos) {
    if (!histos) return nullptr;
    auto new_histo = histos->insert({name, {.0, config.beta(), 100}});
    return &(new_histo.first->second);
  }
*/

  move_worm::move_worm(int block_index, int block_size, std::string const &block_name, qmc_data &data, wl_data& data_wl, bool yes_worm, mc_tools::random_generator &rng, histo_map_t *histos)
  : data_wl(data_wl),
    rng(rng),
    worm_insert(block_index, block_size, block_name, data, data_wl, yes_worm, rng, histos),
    worm_remove(block_index, block_size, block_name, data, data_wl, yes_worm, rng, histos)
  {}
  mc_weight_t move_worm::attempt() {
 

//    std::cout<<"Attempting a worm move !!!, is worm container empty?"<<data_wl.no_worm()<<std::endl;
    if(data_wl.no_worm())
    {

    std::cout<<"Attempting a worm insertion !!!"<<std::endl;
    return worm_insert.attempt();

    }
    else
    {
    std::cout<<"Attempting a worm remove !!!"<<std::endl;
    
    return worm_remove.attempt();

    }
  }

  mc_weight_t move_worm::accept() {

    if(data_wl.no_worm())
    {

    return worm_insert.accept();

    }
    else
    {
    
    return worm_remove.accept();

    }
  }

  void move_worm::reject() {
  
    if(data_wl.no_worm())
    {

    worm_insert.reject();

    }
    else
    {
    
    worm_remove.reject();

    }
  }
}
