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

#pragma once

#include "../common.hpp"
#include "../interfaces/gsl_integration_wrap.hpp"

#include <triqs/gfs.hpp>
#include <triqs/utility/first_include.hpp>
#include <triqs/utility/variant.hpp>
#include <vector>

using namespace triqs::gfs;

namespace keldy::impurity_oneband {

class CPP2PY_IGNORE model_param_t {
 public:
  double beta = 1.0;
  double bias_V = 0.0;
  double eps_d = 0.0;
  double Gamma = 1.0;
  double time_max = +100.0; // (time_limit_min = - time_limit_max)
  int nr_time_points_gf = 1000;
  double alpha = 0.0;
  std::string bath_type = "flatband_fft";
};

// fake function to get cpp2py to create adaptor for model_param_t
CPP2PY_ARG_AS_DICT inline void fake([[maybe_unused]] model_param_t const &temp){};

// ***************************

/// Index of the Keldysh Green Function
class gf_index_t {
 public:
  contour_pt_t contour{};
  spin_t spin = up;
  orbital_t orbital = 0;

  /// Constructor: (time, spin, keldysh_idx)
  //TODO: how to deal with optionnal parameters? (Corentin)
  gf_index_t() = default;
  gf_index_t(time_real_t time_, int spin_, int k_idx_) : contour{time_, keldysh_idx_t(k_idx_)}, spin(spin_t(spin_)) {}
  gf_index_t(time_real_t time_, int spin_, int k_idx_, int timesplit_n_)
     : contour{time_, keldysh_idx_t(k_idx_), timesplit_n_}, spin(spin_t(spin_)) {}
  gf_index_t(time_real_t time_, int spin_, int k_idx_, int timesplit_n_, int orbital_)
     : contour{time_, keldysh_idx_t(k_idx_), timesplit_n_}, spin(spin_t(spin_)), orbital(orbital_t(orbital_)) {}
};

/// Define a total ordering over gf_index_t values. Useful for sorting.
bool operator<(gf_index_t const &a, gf_index_t const &b);

// no time-split comparison in equality
inline bool equivalent_without_timesplit(const gf_index_t &lhs, const gf_index_t &rhs) {
  return (lhs.contour.time == rhs.contour.time) && (lhs.contour.k_idx == rhs.contour.k_idx) && (lhs.spin == rhs.spin)
     && (lhs.orbital == rhs.orbital);
}

inline std::ostream &operator<<(std::ostream &os, gf_index_t const &m) {
  return os << "{" << m.contour.time << ", " << m.contour.k_idx << ", " << m.spin << ", " << m.orbital << "}";
}

// ***************************

/// Defines model throuh non-interacting Green function g_lesser / g_greater
class g0_model {
 public:
  /// Lesser Green function $G^{<}_{\sigma}(t)$; block spin $\sigma$ {up, down}
  block_gf<retime, matrix_valued> g0_lesser;
  array<double, 2> lesser_ft_error = {{0, 0}, {0, 0}};

  /// Greater Green function $G^{>}_{\sigma}(t)$; block spin $\sigma$ {up, down}
  block_gf<retime, matrix_valued> g0_greater;
  array<double, 2> greater_ft_error = {{0, 0}, {0, 0}};

  explicit g0_model(model_param_t const &parameters, bool with_leads);

  /// make dot g0
  void make_g0_by_fft();
  void make_g0_by_contour(double left_turn_pt, double right_turn_pt);
  void make_flat_band_analytic();

  model_param_t param_; // g0_keldysh_contour_t will need access to alpha
  bool const contain_leads;

  dcomplex g0_R_dot(dcomplex omega) {
    return 1.0 / (omega - param_.eps_d - bath_hybrid_R_left(omega) - bath_hybrid_R_right(omega));
  }

  dcomplex g0_A_dot(dcomplex omega) {
    return 1.0 / (omega - param_.eps_d - bath_hybrid_A_left(omega) - bath_hybrid_A_right(omega));
  }

  dcomplex bath_hybrid_R_left(dcomplex omega) { return bath_hybrid_R_left_(omega); }
  dcomplex bath_hybrid_A_left(dcomplex omega) { return std::conj(bath_hybrid_R_left_(std::conj(omega))); }
  dcomplex bath_hybrid_K_left(dcomplex omega) {
    double mu_left = -param_.bias_V / 2;
    return -(2 * n_fermi(omega - mu_left, param_.beta) - 1.) * (bath_hybrid_R_left(omega) - bath_hybrid_A_left(omega));
  }

  dcomplex bath_hybrid_R_right(dcomplex omega) { return bath_hybrid_R_right_(omega); }
  dcomplex bath_hybrid_A_right(dcomplex omega) { return std::conj(bath_hybrid_R_right_(std::conj(omega))); }
  dcomplex bath_hybrid_K_right(dcomplex omega) {
    double mu_right = +param_.bias_V / 2;
    return -(2 * n_fermi(omega - mu_right, param_.beta) - 1.)
       * (bath_hybrid_R_right(omega) - bath_hybrid_A_right(omega));
  };

 private:
  std::function<dcomplex(dcomplex)> bath_hybrid_R_left_ = []([[maybe_unused]] dcomplex omega) { return 0; };
  std::function<dcomplex(dcomplex)> bath_hybrid_R_right_ = []([[maybe_unused]] dcomplex omega) { return 0; };
};

// ***************************

/// Adapt g0_lesser and g0_greater into Green function on Keldysh contour
struct g0_keldysh_contour_t {
  g0_model model;
  g0_keldysh_contour_t(g0_model model_) : model(std::move(model_)){};

  /// return $g^{ab}(t,t')$ in contour basis from $g^<, g^>$ functions
  dcomplex operator()(gf_index_t const &a, gf_index_t const &b, bool internal_point = true) const {
    if (a.spin != b.spin) {
      return 0.0; //  g0 is diagonal in spin
    }
    int time_order = compare_3way(a.contour, b.contour); // Orders on Keldysh Contour

    if (time_order > 0) {
      return model.g0_greater[a.spin](a.contour.time - b.contour.time)(a.orbital, b.orbital);
    }
    if (time_order < 0) {
      return model.g0_lesser[a.spin](a.contour.time - b.contour.time)(a.orbital, b.orbital);
    }
    // if at equal contour points (incl. time-split)
    if (a.orbital != b.orbital) {
      return model.g0_lesser[a.spin](0.0)(a.orbital, b.orbital);
    }
    // if at equal contour points (incl. time-split) AND equal orbitals
    return model.g0_lesser[a.spin](0.0)(a.orbital, a.orbital)
       - static_cast<int>(internal_point) * 1_j * model.param_.alpha;
  }
};

} // namespace keldy::impurity_oneband
