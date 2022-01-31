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

#include "./G_tau_worm.hpp"

namespace triqs_cthyb {

  using namespace triqs::gfs;
  using namespace triqs::mesh;

  measure_G_tau_worm::measure_G_tau_worm(qmc_data const &data, wl_data& data_wl,int n_tau, gf_struct_t const &gf_struct, container_set_t &results)
     : data(data), data_wl(data_wl), average_sign(0) {
    results.G_tau_accum = block_gf<imtime, G_target_t>({data.config.beta(), Fermion, n_tau}, gf_struct);
    G_tau.rebind(*results.G_tau_accum);
    G_tau() = 0.0;

    results.asymmetry_G_tau = block_gf{G_tau};
    asymmetry_G_tau.rebind(*results.asymmetry_G_tau);
    average_visit=0;
  }

  void measure_G_tau_worm::accumulate(mc_weight_t s) {
    s *= data.atomic_reweighting;

  /*  
    for (auto block_idx : range(G_tau.size())) {
      foreach (data.dets[block_idx], [this, s, block_idx](op_t const &x, op_t const &y, det_scalar_t M) {
        // beta-periodicity is implicit in the argument, just fix the sign properly
        auto val    = (y.first >= x.first ? s : -s) * M;
        double dtau = double(y.first - x.first);
        this->G_tau[block_idx][closest_mesh_pt(dtau)](y.second, x.second) += val;
      })
        ;
    }
*/
    if(!data_wl.no_worm())
    {
    //average_sign += s;
    //average_visit += 1;
        auto block_idx = data_wl.get_block_index(0); 
        // beta-periodicity is implicit in the argument, just fix the sign properly
        auto val    = (data_wl.get_time_dag(0) >= data_wl.get_time(0) ? s : -s);
        double dtau = double(data_wl.get_time_dag(0) - data_wl.get_time(0) );
        G_tau[block_idx][closest_mesh_pt(dtau)] += val;
    }
    else
    {
    average_sign += s;
    average_visit += 1;
    }
  }

  void measure_G_tau_worm::collect_results(mpi::communicator const &c) {

    G_tau        = mpi::all_reduce(G_tau, c);
    average_sign = mpi::all_reduce(average_sign, c);

    for (auto &G_tau_block : G_tau) {
      double beta = G_tau_block.mesh().domain().beta;
      G_tau_block /= -real(average_sign) * beta * G_tau_block.mesh().delta();

      // Multiply first and last bins by 2 to account for full bins
      int last = G_tau_block.mesh().size() - 1;
      G_tau_block[0] *= 2;
      G_tau_block[last] *= 2;

      // Set 1/iw behaviour of tails in G_tau to avoid problems when taking FTs later
      auto d = max_element(abs(G_tau_block[0] + G_tau_block[last] + 1));
      if (d > 1e-2 && c.rank() == 0)
        std::cerr << "WARNING: Tau discontinuity of G_tau deviates appreciably from -1\n     .... max_element |g(0) + g(beta) + 1| = " << d << "\n";

      G_tau_block[last] = -1. - G_tau_block[0]; // Enforce 1/iw discontinuity (nb. matrix eq.)
    }

    std::cout<<"Sign Z in worm meas: "<<average_sign<<" and num vistt Z : "<<average_visit<<" with average: "<<double(average_sign/double(average_visit))<<std::endl;

    // We enforce the fundamental Green function property G(tau)[i,j] = G(tau)*[j,i]
    // and store the symmetry violation separately
    asymmetry_G_tau = make_hermitian(G_tau) - G_tau;
    G_tau           = G_tau + asymmetry_G_tau;
  }

} // namespace triqs_cthyb
