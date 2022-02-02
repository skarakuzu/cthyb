/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, H. U.R. Strand, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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
#pragma once

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

#include "../qmc_data.hpp"
#include "../container_set.hpp"

#include "../wl_data.hpp"

namespace triqs_cthyb {

  using namespace triqs::gfs;
  using namespace triqs::mesh;

  // Measure imaginary time Green's function (all blocks)
  class measure_G_tau_worm {

    public:
    measure_G_tau_worm(qmc_data const &data, wl_data& data_wl, int n_tau, gf_struct_t const &gf_struct, container_set_t &results);
    void accumulate(mc_weight_t s);
    void collect_results(mpi::communicator const &c);

    private:
    qmc_data const &data;
    wl_data  &data_wl;
    mc_weight_t average_sign_Z;
    mc_weight_t average_sign_Gup;
    mc_weight_t average_sign_Gdwn;
    G_tau_G_target_t::view_type G_tau;
    G_tau_G_target_t::view_type asymmetry_G_tau;
    int average_visit_Z;
    int average_visit_Gup;
    int average_visit_Gdwn;
  };

} // namespace triqs_cthyb
