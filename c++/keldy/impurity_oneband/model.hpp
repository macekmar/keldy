/*******************************************************************************
 *
 * keldy
 *
 * Copyright (c) 2019-2020 The Simons Foundation
 *   authors: Philipp Dumitrescu, Corentin Bertrand
 *
 * keldy is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * keldy is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * keldy. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include "../common.hpp"
#include "../interfaces/gsl_integration_wrap.hpp"

#include <triqs/gfs.hpp>
#include <triqs/utility/first_include.hpp>
#include <variant>
#include <vector>

namespace keldy::impurity_oneband {

using namespace triqs::gfs;

class CPP2PY_IGNORE model_param_t {
 public:
  double beta = 1.0;
  double bias_V = 0.0;
  double eps_d = 0.0;
  double Gamma = 1.0;
  double alpha = 0.0;
  double half_bandwidth = 2.;
  std::string bath_type = "semicircle";
  // ********
  double time_max = +100.0; // (time_limit_min = - time_limit_max)
  int nr_time_points_gf = 1000;
  std::string ft_method = "fft"; // or "contour" or "analytic"
};

// fake function to get cpp2py to create adaptor for model_param_t
CPP2PY_ARG_AS_DICT inline void fake([[maybe_unused]] model_param_t const &temp){}

void h5_write(h5::group &h5group, std::string const &subgroup_name, model_param_t const &c);
void h5_read(h5::group &h5group, std::string const &subgroup_name, model_param_t &c);

// *****************************************************************************

/// Index of the Keldysh Green Function
class gf_index_t {
 public:
  contour_pt_t contour{};
  spin_t spin = up;
  orbital_t orbital = 0;

  /// Constructor: (time, spin, keldysh_idx)
  // TODO: how to deal with optionnal parameters? (Corentin)
  // TODO: Clean this up (only need cases where enums not construcitble from python)
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

// *****************************************************************************

/// Defines model throuh non-interacting Green function g_lesser / g_greater
class g0_model_omega {
 private:
  model_param_t param_;
  std::function<dcomplex(dcomplex)> bath_hybrid_R_left_ = []([[maybe_unused]] dcomplex omega) { return 0; };
  std::function<dcomplex(dcomplex)> bath_hybrid_R_right_ = []([[maybe_unused]] dcomplex omega) { return 0; };

  dcomplex g0_dot(dcomplex omega, bool is_lesser) const {
    auto bath_hybrid_left = bath_hybrid_R_left(omega) - bath_hybrid_A_left(omega);
    auto bath_hybrid_right = bath_hybrid_R_right(omega) - bath_hybrid_A_right(omega);
    int sign = 2 * static_cast<int>(is_lesser) - 1; // +1 if lesser, -1 if greater
    return -sign * g0_dot_R(omega)
       * (n_fermi(sign * (omega - mu_left()), param_.beta) * bath_hybrid_left
          + n_fermi(sign * (omega - mu_right()), param_.beta) * bath_hybrid_right)
       * g0_dot_A(omega);
  }

  dcomplex g0_rightlead_dot(dcomplex omega, bool is_lesser) const {
    auto bath_hybrid_right = bath_hybrid_R_right(omega) - bath_hybrid_A_right(omega);
    int sign = 2 * static_cast<int>(is_lesser) - 1; // +1 if lesser, -1 if greater
    return -sign * n_fermi(sign * (omega - mu_right()), param_.beta) * bath_hybrid_right * g0_dot_A(omega)
       + bath_hybrid_R_right(omega) * g0_dot(omega, is_lesser);
  }

 public:
  g0_model_omega() = default;
  explicit g0_model_omega(model_param_t const &parameters);

  model_param_t get_param() const { return param_; }

  double get_param_alpha() const { return param_.alpha; }

  double mu_left() const { return -param_.bias_V / 2; }
  double mu_right() const { return +param_.bias_V / 2; }

  dcomplex g0_dot_R(dcomplex omega) const {
    return 1.0 / (omega - param_.eps_d - bath_hybrid_R_left(omega) - bath_hybrid_R_right(omega));
  }
  dcomplex g0_dot_A(dcomplex omega) const {
    return 1.0 / (omega - param_.eps_d - bath_hybrid_A_left(omega) - bath_hybrid_A_right(omega));
  }
  dcomplex g0_dot_lesser(dcomplex omega) const { return g0_dot(omega, true); }
  dcomplex g0_dot_greater(dcomplex omega) const { return g0_dot(omega, false); }
  dcomplex g0_dot_keldysh(keldysh_idx_t a, keldysh_idx_t b, dcomplex omega) const {
    if (a < b) {
      return g0_dot(omega, true); // lesser
    }
    if (a > b) {
      return g0_dot(omega, false); // greater
    }
    // a == b
    if (a == forward) {
      return g0_dot(omega, true) + g0_dot_R(omega); // time-ordered
    }
    // a == backward
    return g0_dot(omega, false) - g0_dot_R(omega); // anti time-ordered
  }

  dcomplex bath_hybrid_R_left(dcomplex omega) const { return bath_hybrid_R_left_(omega); }
  dcomplex bath_hybrid_A_left(dcomplex omega) const { return std::conj(bath_hybrid_R_left_(std::conj(omega))); }
  dcomplex bath_hybrid_K_left(dcomplex omega) const {
    return -(2 * n_fermi(omega - mu_left(), param_.beta) - 1.)
       * (bath_hybrid_R_left(omega) - bath_hybrid_A_left(omega));
  }

  dcomplex bath_hybrid_R_right(dcomplex omega) const { return bath_hybrid_R_right_(omega); }
  dcomplex bath_hybrid_A_right(dcomplex omega) const { return std::conj(bath_hybrid_R_right_(std::conj(omega))); }
  dcomplex bath_hybrid_K_right(dcomplex omega) const {
    return -(2 * n_fermi(omega - mu_right(), param_.beta) - 1.)
       * (bath_hybrid_R_right(omega) - bath_hybrid_A_right(omega));
  }

  dcomplex g0_rightlead_dot_lesser(dcomplex omega) const { return g0_rightlead_dot(omega, true); }
  dcomplex g0_rightlead_dot_greater(dcomplex omega) const { return g0_rightlead_dot(omega, false); }

  friend void h5_write(h5::group &h5group, std::string const &subgroup_name, g0_model_omega const &c);
  friend void h5_read(h5::group &h5group, std::string const &subgroup_name, g0_model_omega &c);
};

// *****************************************************************************

class g0_model {
 public:
  g0_model_omega model_omega{};
  bool make_dot_lead = false;

  /// Lesser Green function $G^{<}_{\sigma}(t)$; block spin $\sigma$ {up, down}
  block_gf<retime, matrix_valued> g0_lesser;
  gf<retime, matrix_valued> lesser_ft_error;

  /// Greater Green function $G^{>}_{\sigma}(t)$; block spin $\sigma$ {up, down}
  block_gf<retime, matrix_valued> g0_greater;
  gf<retime, matrix_valued> greater_ft_error;

  g0_model() = default;
  g0_model(g0_model_omega model_omega_, bool make_dot_lead_);

  static std::string hdf5_format() { return "KELDY_G0Model"; }

  friend void h5_write(h5::group &h5group, std::string const &subgroup_name, g0_model const &c) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write_attribute(grp, "Format", g0_model::hdf5_format());
    h5_write(grp, "model_omega", c.model_omega);
    h5_write(grp, "make_dot_lead", c.make_dot_lead);
    h5_write(grp, "g0_lesser", c.g0_lesser);
    h5_write(grp, "lesser_ft_error", c.lesser_ft_error);
    h5_write(grp, "g0_greater", c.g0_greater);
    h5_write(grp, "greater_ft_error", c.greater_ft_error);
  }

  // Function that read all containers to hdf5 file
  friend void h5_read(h5::group &h5group, std::string const &subgroup_name, g0_model &c) {
    auto grp = h5group.open_group(subgroup_name);
    h5_read(grp, "model_omega", c.model_omega);
    h5_read(grp, "make_dot_lead", c.make_dot_lead);
    h5_read(grp, "g0_lesser", c.g0_lesser);
    h5_read(grp, "lesser_ft_error", c.lesser_ft_error);
    h5_read(grp, "g0_greater", c.g0_greater);
    h5_read(grp, "greater_ft_error", c.greater_ft_error);
  }

  // Function that read all containers to hdf5 file
  CPP2PY_IGNORE
  static g0_model h5_read_construct(h5::group h5group, std::string const &subgroup_name) {
    g0_model c{};
    h5_read(h5group, subgroup_name, c);
    return c;
  }

 private:
  /// make dot g0
  void make_g0_by_fft();
  void make_g0_by_finite_contour(std::vector<double> pts);
  void make_g0_by_contour(double left_turn_pt, double right_turn_pt);
  void make_flat_band_analytic();
};

// *****************************************************************************

/// Adapt g0_lesser and g0_greater into Green function on Keldysh contour
// TODO: why don't we make this a method of model?
struct g0_keldysh_contour_t {

  g0_model model;
  g0_keldysh_contour_t(g0_model model_) : model(std::move(model_)){};

  /// return $g^{ab}(t,t')$ in contour basis from $g^<, g^>$ functions
  dcomplex operator()(gf_index_t const &a, gf_index_t const &b, bool internal_point = true) const {
    using namespace std::complex_literals;

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
       - static_cast<int>(internal_point) * 1.0i * model.model_omega.get_param_alpha();
  }

  double get_time_max() const { return model.g0_lesser[up].mesh().x_max(); };
};

} // namespace keldy::impurity_oneband
