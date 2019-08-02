//******************************************************************************
//
// keldy
//
// Copyright (C) 2019, The Simons Foundation
// authors: Philipp Dumitrescu
//
// keldy is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// keldy is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// keldy. If not, see <http://www.gnu.org/licenses/>.
//
//******************************************************************************

#include "model.hpp"
#include <triqs/gfs.hpp>

namespace keldy::impurity_oneband {

g0_model::g0_model(model_param_t const &parameters) : param_(parameters) { 
  using namespace triqs::gfs;

  auto time_mesh = gf_mesh<retime>({-param_.time_max, param_.time_max, param_.nr_time_points_gf});
  auto freq_mesh = make_adjoint_mesh(time_mesh);

  gf<refreq, scalar_valued> g0_lesser_omega{freq_mesh};
  gf<refreq, scalar_valued> g0_greater_omega{freq_mesh};

  // Define local Fermi funciton; reads in beta
  auto nFermi = [&](double omega) {
    if (omega > 0){
      auto y = std::exp(- param_.beta * omega);
      return y / (1. + y);
    }
    return  1.0 / (std::exp(param_.beta * omega) + 1);
  };

  for(auto w : freq_mesh){
    g0_lesser_omega[w] = 1_j * param_.Gamma * (nFermi(w + param_.bias_V / 2) + nFermi(w - param_.bias_V / 2))
          / ((w - param_.eps_d) * (w - param_.eps_d) + std::pow(param_.Gamma, 2));
    // g0_lesser_omega()[down] = g0_lesser_omega()[up];

    g0_greater_omega[w] = 1_j * param_.Gamma * (nFermi(w + param_.bias_V / 2) + nFermi(w - param_.bias_V / 2) - 2)
          / ((w - param_.eps_d) * (w - param_.eps_d) + std::pow(param_.Gamma, 2));
    // g0_greater_omega(w)[down] = g0_lesser_omega(w)[up];

  }

  gf<retime, scalar_valued> g0_lesser_up = make_gf_from_fourier(g0_lesser_omega, time_mesh);
  gf<retime, scalar_valued> g0_greater_up = make_gf_from_fourier(g0_greater_omega, time_mesh);

  // Since Spin up and down are currently identical
  g0_lesser = make_block_gf<retime, scalar_valued>({"up", "down"}, {g0_lesser_up, g0_lesser_up});
  g0_greater = make_block_gf<retime, scalar_valued>({"up", "down"}, {g0_greater_up, g0_greater_up});
}




dcomplex g0_keldysh_contour_t::operator()(gf_index_t const &a, gf_index_t const &b, bool use_lesser_at_eq_points) const {
  if (a.spin != b.spin) {
    return 0.0; //  g0 is diagonal in spin
  }

  // At equal contour points, we need to point split:
  // For self-contraction: use g0_lesser (~ c^\dag c), since V is normal ordered
  // Multiple vertices at equal time: need to correctly swap relative g0_greater / g_lesser
  // (Float is ok in == since we are not doing any operations on times),
  if (a.time == b.time && a.k_idx == b.k_idx) {
    if (!use_lesser_at_eq_points){
      return model.g0_greater[a.spin](0.0);
    }
    return model.g0_lesser[a.spin](0.0);
  }

  // mapping to g^< / g^> [0 = forward contour; 1 = backward contour]
  //  a    b    (a.time > b.time)   L/G ?
  //  0    0           1             G
  //  0    0           0             L
  //  1    1           1             L
  //  1    1           0             G
  //
  //  0    1           *             L
  //  1    0           *             G

  bool is_greater = (a.k_idx == b.k_idx ? (a.k_idx xor (a.time > b.time)) : a.k_idx);

  if (is_greater) {
    return model.g0_greater[a.spin](a.time - b.time);
  }
  return model.g0_lesser[a.spin](a.time - b.time);
}

} // namespace keldy::impurity_oneband
