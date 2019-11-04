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

bool operator<(gf_index_t const &a, gf_index_t const &b) {
  if (a.k_idx != b.k_idx) {
    return (a.k_idx < b.k_idx);
  } else {
    if (a.time != b.time) {
      return (a.k_idx == forward) ? a.time < b.time : -a.time < -b.time;
    } else {
      if (a.timesplit_n != b.timesplit_n) {
        return (a.k_idx == forward) ? a.timesplit_n < b.timesplit_n : -a.timesplit_n < -b.timesplit_n;
      } else {
        return a.spin < b.spin;
      }
    }
  }
}

g0_model::g0_model(model_param_t const &parameters) : param_(parameters) {
  if (param_.bath_type == "semicircle") {
    make_semicircular_model();
  } else if (param_.bath_type == "flatband") {
    make_flat_band();
  } else {
    TRIQS_RUNTIME_ERROR << "bath_type not defined";
  }
}

void g0_model::make_semicircular_model() {
  using namespace triqs::gfs;

  auto time_mesh = gf_mesh<retime>({-param_.time_max, param_.time_max, param_.nr_time_points_gf});
  auto freq_mesh = make_adjoint_mesh(time_mesh);

  gf<refreq, scalar_valued> g0_lesser_omega{freq_mesh};
  gf<refreq, scalar_valued> g0_greater_omega{freq_mesh};

  // Define local Fermi funciton; reads in beta
  auto nFermi = [&](double omega) {
    if (omega > 0) {
      auto y = std::exp(-param_.beta * omega);
      return y / (1. + y);
    }
    return 1.0 / (std::exp(param_.beta * omega) + 1);
  };

  // Retarded self energy with semi circular sigma dos (linear chain).
  auto sigma_linear_chain = [](double omega) -> dcomplex {
    omega = omega / 2.;
    if (std::abs(omega) < 1) {
      return dcomplex{omega, -std::sqrt(1 - omega * omega)};
    }
    if (omega > 1) {
      return omega - std::sqrt(omega * omega - 1);
    }
    return omega + std::sqrt(omega * omega - 1);
  };

  // The non interacting dot GF's in frequency (2*2 matrix with Keldysh indices)
  // From XW computation
  auto G0_dd_w = [&](double omega) {
    dcomplex gr = 1 / (omega - param_.eps_d - 2 * param_.Gamma * sigma_linear_chain(omega));
    dcomplex fac = 2_j * gr * conj(gr);
    dcomplex gam_gg = -1 * param_.Gamma * imag(sigma_linear_chain(omega)) * fac;
    dcomplex temp = (nFermi(omega - param_.bias_V / 2) + nFermi(omega - param_.bias_V / 2)) * gam_gg;
    dcomplex temp2 = temp - 2 * gam_gg;
    return array<dcomplex, 2>{{gr + temp, temp}, {temp2, -conj(gr) + temp}};
  };

  for (auto w : freq_mesh) {
    auto g0_dd = G0_dd_w(w);

    g0_lesser_omega[w] = g0_dd(0, 1);
    g0_greater_omega[w] = g0_dd(1, 0);
  }

  gf<retime, scalar_valued> g0_lesser_up = make_gf_from_fourier(g0_lesser_omega, time_mesh);
  gf<retime, scalar_valued> g0_greater_up = make_gf_from_fourier(g0_greater_omega, time_mesh);

  // Since Spin up and down are currently identical
  g0_lesser = make_block_gf<retime, scalar_valued>({"up", "down"}, {g0_lesser_up, g0_lesser_up});
  g0_greater = make_block_gf<retime, scalar_valued>({"up", "down"}, {g0_greater_up, g0_greater_up});
}

void g0_model::make_flat_band() {
  using namespace triqs::gfs;

  auto time_mesh = gf_mesh<retime>({-param_.time_max, param_.time_max, param_.nr_time_points_gf});
  auto freq_mesh = make_adjoint_mesh(time_mesh);

  gf<refreq, scalar_valued> g0_lesser_omega{freq_mesh};
  gf<refreq, scalar_valued> g0_greater_omega{freq_mesh};

  // Define local Fermi funciton; reads in beta
  auto nFermi = [&](double omega) {
    if (omega > 0) {
      auto y = std::exp(-param_.beta * omega);
      return y / (1. + y);
    }
    return 1.0 / (std::exp(param_.beta * omega) + 1);
  };

  auto G0_dd_w = [&](double w) {
    double we = w - param_.eps_d;
    auto R = 0.5 / (we + 1_j * param_.Gamma);
    auto A = 0.5 / (we - 1_j * param_.Gamma);
    auto K = 1_j * param_.Gamma / (we * we + param_.Gamma * param_.Gamma)
       * (nFermi(w - param_.bias_V / 2) + nFermi(w + param_.bias_V / 2) - 1.);
    return array<dcomplex, 2>{{K + R + A, K - R + A}, {K + R - A, K - R - A}};
  };

  for (auto w : freq_mesh) {
    auto g0_dd = G0_dd_w(w);

    g0_lesser_omega[w] = g0_dd(0, 1);
    g0_greater_omega[w] = g0_dd(1, 0);

    // g0_lesser_omega[w] = 1_j * param_.Gamma * (nFermi(w + param_.bias_V / 2) + nFermi(w - param_.bias_V / 2))
    //    / ((w - param_.eps_d) * (w - param_.eps_d) + std::pow(param_.Gamma, 2));
    // // g0_lesser_omega()[down] = g0_lesser_omega()[up];

    // g0_greater_omega[w] = 1_j * param_.Gamma * (nFermi(w + param_.bias_V / 2) + nFermi(w - param_.bias_V / 2) - 2)
    //    / ((w - param_.eps_d) * (w - param_.eps_d) + std::pow(param_.Gamma, 2));
    // // g0_greater_omega(w)[down] = g0_lesser_omega(w)[up];
  }

  gf<retime, scalar_valued> g0_lesser_up = make_gf_from_fourier(g0_lesser_omega, time_mesh);
  gf<retime, scalar_valued> g0_greater_up = make_gf_from_fourier(g0_greater_omega, time_mesh);

  // Since Spin up and down are currently identical
  g0_lesser = make_block_gf<retime, scalar_valued>({"up", "down"}, {g0_lesser_up, g0_lesser_up});
  g0_greater = make_block_gf<retime, scalar_valued>({"up", "down"}, {g0_greater_up, g0_greater_up});
}

// *****

/// Adaptor: return $g^{ab}(t,t')$ in contour basis from $g^<, g^>$ functions
/// Mapping takes care of correct keldysh ordering [0 = forward contour; 1 = backward contour]
///  a    b    (a.time > b.time)   L/G ?
///  0    0           1             G
///  0    0           0             L
///  1    1           1             L
///  1    1           0             G
///
///  0    1           *             L
///  1    0           *             G
///
/// At equal times and Keldysh index we need to point-split times according to external integer (time-ordering).
/// For self contractions (no external ordering), use $g^< \sim c^\dag c$ since this is normal ordering defined by V
dcomplex g0_keldysh_contour_t::operator()(gf_index_t const &a, gf_index_t const &b, bool internal_point) const {
  if (a.spin != b.spin) {
    return 0.0; //  g0 is diagonal in spin
  }

  // Takes care of 01 & 10 case
  bool is_greater = a.k_idx;

  if (a.k_idx == b.k_idx) {
    if (a.time != b.time) { // generic case (table above)
      is_greater = a.k_idx xor (a.time > b.time);
    } else {
      // Need to time split from external ordering unless it is a self-contraction (tadpole)
      if (a.timesplit_n != b.timesplit_n) {
        is_greater = a.k_idx xor (a.timesplit_n > b.timesplit_n);
      } else {
        // if internal index use alpha. Else use from external
        return model.g0_lesser[a.spin](0.0) - static_cast<long>(internal_point) * 1_j * model.param_.alpha;
      }
    }
  }

  if (is_greater) {
    return model.g0_greater[a.spin](a.time - b.time);
  }
  return model.g0_lesser[a.spin](a.time - b.time);
}

} // namespace keldy::impurity_oneband
