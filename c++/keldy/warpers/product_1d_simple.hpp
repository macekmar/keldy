/*******************************************************************************
 *
 * keldy
 *
 * Copyright (c) 2020 The Simons Foundation
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
#include "warpers_common.hpp"
#include <algorithm>
#include <any>
#include <functional>
#include <iomanip>
#include <itertools/omp_chunk.hpp>
#include <numeric>
#include <triqs/gfs.hpp>
#include <memory>
#include <triqs/utility/exceptions.hpp>
#include <triqs/utility/macros.hpp>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

namespace keldy::warpers {

using gf_t = triqs::gfs::gf<triqs::gfs::retime, triqs::gfs::scalar_real_valued>;

class warper_product_1d_simple_t {
 private:
  double f1_integrated_normalization = 1.0;

  std::function<double(double)> f1 = [](double /*u*/) { return 1.0; };
  std::function<double(double)> f1_integrated = [](double u) { return u; };
  std::function<double(double)> f1_integrated_inverse = [](double l) { return l; };

  constexpr static double domain_u_min = 0.0; // always start at 0.0
  double domain_u_max = 1.0;                  // "t_max"

  constexpr static double codomain_l_min = 0.0;
  constexpr static double codomain_l_max = 1.0;

 public:
  warper_product_1d_simple_t() = default;

  // Constructor to use if f1, f1_integrated and f1_integrated_inverse can be provided analytically
  warper_product_1d_simple_t(std::function<double(double)> f1_, std::function<double(double)> f1_integrated_,
                             std::function<double(double)> f1_integrated_inverse_, double domain_u_max_,
                             bool do_domain_checks = true)
     : f1(std::move(f1_)),
       f1_integrated(std::move(f1_integrated_)),
       f1_integrated_inverse(std::move(f1_integrated_inverse_)),
       domain_u_max(domain_u_max_) {
    // Basic Check:
    if (!(domain_u_min < domain_u_max)) {
      TRIQS_RUNTIME_ERROR << "coord_domain_u coordinates must be ordered.";
    }

    if (do_domain_checks) {
      double eps = 1e-12;
      bool domain_match = (std::abs(f1_integrated(domain_u_min) - codomain_l_min) < eps)
         && (std::abs(f1_integrated(domain_u_max) - codomain_l_max) < eps)
         && (std::abs(f1_integrated_inverse(codomain_l_min) - domain_u_min) < eps)
         && (std::abs(f1_integrated_inverse(codomain_l_max) - domain_u_max) < eps);

      if (!domain_match) {
        TRIQS_RUNTIME_ERROR << "f1_integrated and f1_integrated_inverse do not map domain boundaries correctly.\n"
                            << "* std::abs(f1_integrated(domain_u_min) - codomain_l_min): "
                            << (std::abs(f1_integrated(domain_u_min) - codomain_l_min)) << "\n"
                            << "* std::abs(f1_integrated(domain_u_max) - codomain_l_max): "
                            << (std::abs(f1_integrated(domain_u_max) - codomain_l_max)) << "\n"
                            << "* std::abs(f1_integrated_inverse(codomain_l_min) - domain_u_min): "
                            << (std::abs(f1_integrated_inverse(codomain_l_min) - domain_u_min)) << "\n"
                            << "* std::abs(f1_integrated_inverse(codomain_l_max) - domain_u_max): "
                            << (std::abs(f1_integrated_inverse(codomain_l_max) - domain_u_max)) << "\n";
      }
      // f1 > 0 : for Jacobian consistancy
      // f1_integrated & f1_integrated_inverse should be monotonous
      // f1_integrated & f1_integrated_inverse are inverses

      // check consistency of functions provided
      for (double l_prime : {0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.}) {
        if (std::abs(f1_integrated(f1_integrated_inverse(l_prime)) - l_prime) > 1e-12) {
          TRIQS_RUNTIME_ERROR
             << "Inconsistent functions: f1_integrated_inverse should be the inverse of f1_integrated ("
             << f1_integrated_(f1_integrated_inverse_(l_prime)) << " != " << l_prime << ")";
        }
      }
    }
  }

  // ****************************************

  [[nodiscard]] std::pair<std::vector<double>, double> map_reverse(std::vector<double> const &li_vec) const {
    auto xi_vec = li_vec;
    double jacobian_r = 1.0;
    for (auto &xi : xi_vec) {
      xi = f1_integrated_inverse(xi);
      jacobian_r *= 1.0 / f1(xi); // f1 defined in ui, so evaluate after map
    }
    return std::make_pair(xi_vec, jacobian_r);
  }

  [[nodiscard]] std::pair<std::vector<double>, double> map_forward(std::vector<double> const &ui_vec) const {
    auto xi_vec = ui_vec;
    double jacobian_f = 1.0;
    for (auto &xi : xi_vec) {
      jacobian_f *= f1(xi); // f1 defined in ui, so evaluate before map
      xi = f1_integrated(xi);
    }
    return std::make_pair(xi_vec, jacobian_f);
  }

  [[nodiscard]] std::vector<double> ui_from_li(std::vector<double> const &li_vec) const {
    std::vector<double> ui_vec = li_vec;
    for (auto &ui : ui_vec) {
      ui = f1_integrated_inverse(ui);
    }
    return ui_vec;
  }

  [[nodiscard]] std::vector<double> li_from_ui(std::vector<double> const &ui_vec) const {
    std::vector<double> li_vec = ui_vec;
    for (auto &li : li_vec) {
      li = f1_integrated(li);
    }
    return li_vec;
  }

  [[nodiscard]] double jacobian_reverse(std::vector<double> const &li_vec) const {
    double result = 1.0;
    for (auto const &li : li_vec) {
      result *= 1.0 / f1(f1_integrated_inverse(li));
    }
    return result;
  }

  [[nodiscard]] double jacobian_forward(std::vector<double> const &ui_vec) const {
    double result = 1.0;
    for (auto const &ui : ui_vec) {
      result *= f1(ui);
    }
    return result;
  }

  [[nodiscard]] double operator()(std::vector<double> const &ui_vec) const {
    double result = 1.0;
    for (auto const &ui : ui_vec) {
      result *= f1(ui) / f1_integrated_normalization;
    }
    return result;
  }
};

// *****************

namespace nda = triqs::arrays;

class warper_product_1d_simple_interp_nearest_t {
 private:
  double f1_integrated_normalization = 1.0;

  constexpr static double domain_u_min = 0.0; // always start at 0.0
  double domain_u_max = 1.0;                  // "t_max"

  // the codomain is [0, 1], but internally this codomain is reversed onto [-1, 0]
  constexpr static double codomain_l_min = 0.0;
  constexpr static double codomain_l_max = 1.0;

  int nr_sample_points = 0;

  nda::vector<double> times_u_pts{};
  nda::vector<double> f1_pts{};
  nda::vector<double> f1_integrated_pts{};

  std::shared_ptr<gsl_spline> sp_f_integrated;
  std::shared_ptr<gsl_spline> sp_f_integrated_inverse;

  std::shared_ptr<gsl_interp_accel> acc_f_integrated;
  std::shared_ptr<gsl_interp_accel> acc_f_integrated_inverse;

  // f is constant interpolated to nearest downward: this piggy-backs on the acc lookup for f_integrated
  double f1(double u) const {
    assert((domain_u_min <= u) && (u <= domain_u_max));
    auto i = gsl_interp_accel_find(acc_f_integrated.get(), times_u_pts.data_start(), nr_sample_points, u);
    if (i == times_u_pts.size() - 1) {
      return f1_pts[i];
    };
    return 0.5 * (f1_pts[i] + f1_pts[i + 1]);
  }

  double f1_integrated(double u) const {
    assert((domain_u_min <= u) && (u <= domain_u_max));
    double l = gsl_spline_eval(sp_f_integrated.get(), u, acc_f_integrated.get());
    return -l; // map [-1, 0] to [0, 1] without accuracy loss
  }

  double f1_integrated_inverse(double l) const {
    assert((codomain_l_min <= l) && (l <= codomain_l_max));
    l = -l; // map [0, 1] to [-1, 0] without accuracy loss
    return gsl_spline_eval(sp_f_integrated_inverse.get(), l, acc_f_integrated_inverse.get());
  }

  void integrate_and_inverse() {
    // integrate in the negative direction from u_max.
    int cut_index = nr_sample_points - 1;
    bool found_cut = false;
    f1_integrated_pts[nr_sample_points - 1] = 0.0; // scale co-ordinates later
    for (int i = nr_sample_points - 2; i >= 0; i--) {
      f1_integrated_pts[i] =
         f1_integrated_pts[i + 1] - 0.5 * (f1_pts[i + 1] + f1_pts[i]) * (times_u_pts[i + 1] - times_u_pts[i]);

      found_cut = found_cut || (f1_integrated_pts[i] != 0);
      if (not found_cut) {
        cut_index = i;
      }
    }

    if (not found_cut) {
      TRIQS_RUNTIME_ERROR << "Model function is flat zero.";
    }

    f1_integrated_normalization = -f1_integrated_pts[0]; // > 0

    f1_integrated_pts() /= f1_integrated_normalization;
    f1_pts /= f1_integrated_normalization;
    f1_integrated_pts[0] = -1; // normalization seems to be approximate sometimes

    assert(f1_integrated_pts[0] == -1);
    assert(f1_integrated_pts[cut_index] == 0);
    assert(f1_integrated_pts[nr_sample_points - 1] == 0);
    assert(times_u_pts[0] == 0);
    assert(times_u_pts[nr_sample_points - 1] == domain_u_max);

    sp_f_integrated_inverse = {gsl_spline_alloc(gsl_interp_linear, cut_index + 1), gsl_spline_free};

    gsl_spline_init(sp_f_integrated.get(), times_u_pts.data_start(), f1_integrated_pts.data_start(), nr_sample_points);
    gsl_spline_init(sp_f_integrated_inverse.get(), f1_integrated_pts.data_start(), times_u_pts.data_start(),
                    cut_index + 1);
  }

 public:
  warper_product_1d_simple_interp_nearest_t()
     : sp_f_integrated{gsl_spline_alloc(gsl_interp_linear, nr_sample_points), gsl_spline_free},
       sp_f_integrated_inverse{gsl_spline_alloc(gsl_interp_linear, nr_sample_points), gsl_spline_free},
       acc_f_integrated{gsl_interp_accel_alloc(), gsl_interp_accel_free},
       acc_f_integrated_inverse{gsl_interp_accel_alloc(), gsl_interp_accel_free} {};

  warper_product_1d_simple_interp_nearest_t(nda::vector<double> times_u_pts_, nda::vector<double> f1_pts_)
     : domain_u_max{times_u_pts_[times_u_pts_.size() - 1]},
       nr_sample_points{static_cast<int>(times_u_pts_.size())},
       times_u_pts(std::move(times_u_pts_)),
       f1_pts(std::move(f1_pts_)),
       f1_integrated_pts(nr_sample_points),
       sp_f_integrated{gsl_spline_alloc(gsl_interp_linear, nr_sample_points), gsl_spline_free},
       sp_f_integrated_inverse{gsl_spline_alloc(gsl_interp_linear, nr_sample_points), gsl_spline_free},
       acc_f_integrated{gsl_interp_accel_alloc(), gsl_interp_accel_free},
       acc_f_integrated_inverse{gsl_interp_accel_alloc(), gsl_interp_accel_free} {

    if (f1_pts.size() != nr_sample_points) {
      TRIQS_RUNTIME_ERROR << "times_u_pts and f1_pts should have the same length.";
    }
    if (!(nr_sample_points >= 3)) {
      TRIQS_RUNTIME_ERROR << "Expect nr_sample_points to be >= 3. Given was " << nr_sample_points;
    }
    if (times_u_pts[0] != 0.) {
      TRIQS_RUNTIME_ERROR << "times_u_pts should start at 0";
    }
    ///TODO check strict increasing times_u_pts

    integrate_and_inverse();
  }

  warper_product_1d_simple_interp_nearest_t(std::function<double(double)> const &f1_, double domain_u_max_,
                                            int nr_sample_points_)
     : domain_u_max{domain_u_max_},
       nr_sample_points{nr_sample_points_},
       times_u_pts(nr_sample_points),
       f1_pts(nr_sample_points),
       f1_integrated_pts(nr_sample_points),
       sp_f_integrated{gsl_spline_alloc(gsl_interp_linear, nr_sample_points), gsl_spline_free},
       sp_f_integrated_inverse{gsl_spline_alloc(gsl_interp_linear, nr_sample_points), gsl_spline_free},
       acc_f_integrated{gsl_interp_accel_alloc(), gsl_interp_accel_free},
       acc_f_integrated_inverse{gsl_interp_accel_alloc(), gsl_interp_accel_free} {

    if (!(nr_sample_points >= 3)) {
      TRIQS_RUNTIME_ERROR << "Expect nr_sample_points to be >= 3. Given was " << nr_sample_points;
    }

    // Just do naive linear interpolation on grid
    //#pragma omp parallel for
    for (int i = 0; i < nr_sample_points; i++) {
      times_u_pts[i] = domain_u_min + static_cast<double>(i) / (nr_sample_points - 1) * (domain_u_max - domain_u_min);
    }
    times_u_pts[nr_sample_points - 1] = domain_u_max; // avoid float operation inaccuracy
    for (int i = 0; i < nr_sample_points; i++) {
      f1_pts[i] = f1_(times_u_pts[i]);
    }

    integrate_and_inverse();
  }

  // ****************************************

  [[nodiscard]] std::pair<std::vector<double>, double> map_reverse(std::vector<double> const &li_vec) const {
    auto xi_vec = li_vec;
    double jacobian_r = 1.0;
    for (auto &xi : xi_vec) {
      xi = f1_integrated_inverse(xi);
      jacobian_r *= 1.0 / f1(xi); // f1 defined in ui, so evaluate after map
    }
    return std::make_pair(xi_vec, jacobian_r);
  }

  [[nodiscard]] std::pair<std::vector<double>, double> map_forward(std::vector<double> const &ui_vec) const {
    auto xi_vec = ui_vec;
    double jacobian_f = 1.0;
    for (auto &xi : xi_vec) {
      jacobian_f *= f1(xi); // f1 defined in ui, so evaluate before map
      xi = f1_integrated(xi);
    }
    return std::make_pair(xi_vec, jacobian_f);
  }

  [[nodiscard]] std::vector<double> ui_from_li(std::vector<double> const &li_vec) const {
    std::vector<double> ui_vec = li_vec;
    for (auto &ui : ui_vec) {
      ui = f1_integrated_inverse(ui);
    }
    return ui_vec;
  }

  [[nodiscard]] std::vector<double> li_from_ui(std::vector<double> const &ui_vec) const {
    std::vector<double> li_vec = ui_vec;
    for (auto &li : li_vec) {
      li = f1_integrated(li);
    }
    return li_vec;
  }

  [[nodiscard]] double jacobian_reverse(std::vector<double> const &li_vec) const {
    double result = 1.0;
    for (auto const &li : li_vec) {
      result *= 1.0 / f1(f1_integrated_inverse(li));
    }
    return result;
  }

  [[nodiscard]] double jacobian_forward(std::vector<double> const &ui_vec) const {
    double result = 1.0;
    for (auto const &ui : ui_vec) {
      result *= f1(ui);
    }
    return result;
  }

  [[nodiscard]] double operator()(std::vector<double> const &ui_vec) const {
    double result = 1.0;
    for (auto const &ui : ui_vec) {
      result *= f1(ui) / f1_integrated_normalization;
    }
    return result;
  }
};

class warper_product_1d_simple_interp_hybrid_t {
 private:
  double f1_integrated_normalization = 1.0;

  constexpr static double domain_u_min = 0.0; // always start at 0.0
  double domain_u_max = 1.0;                  // "t_max"

  // the codomain is [0, 1], but internally this codomain is reversed onto [-1, 0]
  constexpr static double codomain_l_min = 0.0;
  constexpr static double codomain_l_max = 1.0;

  int nr_sample_points = 0;

  nda::vector<double> times_u_pts{};
  nda::vector<double> f1_integrated_pts{};

  std::shared_ptr<gsl_spline> sp_f_integrated;
  std::shared_ptr<gsl_spline> sp_f_integrated_inverse;

  std::shared_ptr<gsl_interp_accel> acc_f_integrated;
  std::shared_ptr<gsl_interp_accel> acc_f_integrated_inverse;

  std::function<double(double)> f1 = [](double u [[maybe_unused]]) { return 1.0; };

  double f1_integrated(double u) const {
    assert((domain_u_min <= u) && (u <= domain_u_max));
    double l = gsl_spline_eval(sp_f_integrated.get(), u, acc_f_integrated.get());
    return -l; // map [-1, 0] to [0, 1] without precision loss
  }

  double f1_integrated_inverse(double l) const {
    assert((codomain_l_min <= l) && (l <= codomain_l_max));
    l = -l; // map [0, 1] to [-1, 0] without precision loss
    return gsl_spline_eval(sp_f_integrated_inverse.get(), l, acc_f_integrated_inverse.get());
  }

 public:
  warper_product_1d_simple_interp_hybrid_t()
     : sp_f_integrated{gsl_spline_alloc(gsl_interp_steffen, nr_sample_points), gsl_spline_free},
       sp_f_integrated_inverse{gsl_spline_alloc(gsl_interp_steffen, nr_sample_points), gsl_spline_free},
       acc_f_integrated{gsl_interp_accel_alloc(), gsl_interp_accel_free},
       acc_f_integrated_inverse{gsl_interp_accel_alloc(), gsl_interp_accel_free} {};

  warper_product_1d_simple_interp_hybrid_t(std::function<double(double)> const &f1_, double domain_u_max_,
                                           int nr_sample_points_)
     : domain_u_max{domain_u_max_},
       nr_sample_points{nr_sample_points_},
       times_u_pts(nr_sample_points),
       f1_integrated_pts(nr_sample_points),
       sp_f_integrated{gsl_spline_alloc(gsl_interp_steffen, nr_sample_points), gsl_spline_free},
       sp_f_integrated_inverse{gsl_spline_alloc(gsl_interp_steffen, nr_sample_points), gsl_spline_free},
       acc_f_integrated{gsl_interp_accel_alloc(), gsl_interp_accel_free},
       acc_f_integrated_inverse{gsl_interp_accel_alloc(), gsl_interp_accel_free} {

    if (!(nr_sample_points >= 3)) {
      TRIQS_RUNTIME_ERROR << "Expect nr_sample_points to be >= 3. Given was " << nr_sample_points;
    }

    // Just do naive linear interpolation on grid
    double delta = (domain_u_max - domain_u_min) / (nr_sample_points - 1);
    //#pragma omp parallel for
    for (int i = 0; i < nr_sample_points; i++) {
      times_u_pts[i] = domain_u_min + delta * i;
    }
    times_u_pts[nr_sample_points - 1] = domain_u_max; // avoid float operation inaccuracy

    int cut_index = nr_sample_points - 1;
    bool found_cut = false;

    // Integrate function using Simpsons rule: Evaluate f1 two times for each output point
    f1_integrated_pts[nr_sample_points - 1] = 0.0;       // scale co-ordinates later
    double f2d = f1_(times_u_pts[nr_sample_points - 1]); // temp -- carry over
    for (int i = nr_sample_points - 2; i >= 0; i--) {
      double f0d = f2d;
      double f1d = f1_(0.5 * (times_u_pts[i + 1] + times_u_pts[i]));
      f2d = f1_(times_u_pts[i]);
      double dl = delta / 6. * (f0d + 4 * f1d + f2d);
      f1_integrated_pts[i] = f1_integrated_pts[i + 1] - dl;

      if (dl < 1e-100) {
        TRIQS_RUNTIME_ERROR << "dl=" << dl << " < 1e-100 is too small for a cubic interpolator in double precision.";
      }

      found_cut = found_cut || (f1_integrated_pts[i] != 0);
      if (not found_cut) {
        cut_index = i;
      }
    }

    if (not found_cut) {
      TRIQS_RUNTIME_ERROR << "Model function is flat zero.";
    }

    f1_integrated_normalization = -f1_integrated_pts[0]; // > 0

    f1_integrated_pts /= f1_integrated_normalization;
    f1_integrated_pts[0] = -1; // normalization seems to be approximate sometimes
    f1 = [norm = f1_integrated_normalization, f1_](double u) { return f1_(u) / norm; };

    assert(f1_integrated_pts[0] == -1);
    assert(f1_integrated_pts[cut_index] == 0);
    assert(f1_integrated_pts[nr_sample_points - 1] == 0);
    assert(times_u_pts[0] == 0);
    assert(times_u_pts[nr_sample_points - 1] == domain_u_max);

    sp_f_integrated_inverse = {gsl_spline_alloc(gsl_interp_steffen, cut_index + 1), gsl_spline_free};

    gsl_spline_init(sp_f_integrated.get(), times_u_pts.data_start(), f1_integrated_pts.data_start(), nr_sample_points);
    gsl_spline_init(sp_f_integrated_inverse.get(), f1_integrated_pts.data_start(), times_u_pts.data_start(),
                    cut_index + 1);
  }

  // ****************************************

  [[nodiscard]] std::pair<std::vector<double>, double> map_reverse(std::vector<double> const &li_vec) const {
    auto xi_vec = li_vec;
    double jacobian_r = 1.0;
    for (auto &xi : xi_vec) {
      xi = f1_integrated_inverse(xi);
      jacobian_r *= 1.0 / f1(xi); // f1 defined in ui, so evaluate after map
    }
    return std::make_pair(xi_vec, jacobian_r);
  }

  [[nodiscard]] std::pair<std::vector<double>, double> map_forward(std::vector<double> const &ui_vec) const {
    auto xi_vec = ui_vec;
    double jacobian_f = 1.0;
    for (auto &xi : xi_vec) {
      jacobian_f *= f1(xi); // f1 defined in ui, so evaluate before map
      xi = f1_integrated(xi);
    }
    return std::make_pair(xi_vec, jacobian_f);
  }

  [[nodiscard]] std::vector<double> ui_from_li(std::vector<double> const &li_vec) const {
    std::vector<double> ui_vec = li_vec;
    for (auto &ui : ui_vec) {
      ui = f1_integrated_inverse(ui);
    }
    return ui_vec;
  }

  [[nodiscard]] std::vector<double> li_from_ui(std::vector<double> const &ui_vec) const {
    std::vector<double> li_vec = ui_vec;
    for (auto &li : li_vec) {
      li = f1_integrated(li);
    }
    return li_vec;
  }

  [[nodiscard]] double jacobian_reverse(std::vector<double> const &li_vec) const {
    double result = 1.0;
    for (auto const &li : li_vec) {
      result *= 1.0 / f1(f1_integrated_inverse(li));
    }
    return result;
  }

  [[nodiscard]] double jacobian_forward(std::vector<double> const &ui_vec) const {
    double result = 1.0;
    for (auto const &ui : ui_vec) {
      result *= f1(ui);
    }
    return result;
  }

  [[nodiscard]] double operator()(std::vector<double> const &ui_vec) const {
    double result = 1.0;
    for (auto const &ui : ui_vec) {
      result *= f1(ui) / f1_integrated_normalization;
    }
    return result;
  }
};

// ***************************************************************************************************************

// Maker Functions:

inline warper_product_1d_simple_t make_product_1d_simple_exponential(double domain_u_max, double w_scale) {
  return {[w_scale, domain_u_max](double t) -> double {
            long double norm = -std::expm1(-static_cast<long double>(domain_u_max) / w_scale);
            return std::exp(-t / w_scale) / (w_scale * norm);
          },
          [w_scale, domain_u_max](double t) -> double {
            long double norm = -std::expm1(-static_cast<long double>(domain_u_max) / w_scale);
            return -std::expm1(-t / w_scale) / norm;
          },
          [w_scale, domain_u_max](double l) -> double {
            long double norm = -std::expm1(-static_cast<long double>(domain_u_max) / w_scale);
            return -w_scale * std::log1p(-l * norm);
          },
          domain_u_max};
}

inline warper_product_1d_simple_t make_product_1d_simple_inverse(double domain_u_max, double w_scale) {
  return {[w_scale, domain_u_max](double t) -> double {
            long double norm = std::log1p(static_cast<long double>(domain_u_max) / w_scale);
            return 1.0 / (norm * (w_scale + t));
          },
          [w_scale, domain_u_max](double t) -> double {
            long double norm = std::log1p(static_cast<long double>(domain_u_max) / w_scale);
            return std::log1p(t / w_scale) / norm;
          },
          [w_scale, domain_u_max](double l) -> double {
            long double norm = std::log1p(static_cast<long double>(domain_u_max) / w_scale);
            return w_scale * std::expm1(l * norm);
          },
          domain_u_max};
}

inline warper_product_1d_simple_t make_product_1d_simple_inverse_square(double domain_u_max, double w_scale) {
  return {[w_scale, domain_u_max](double t) -> double {
            long double norm = static_cast<long double>(domain_u_max) / (w_scale + domain_u_max);
            return w_scale / ((w_scale + t) * (w_scale + t) * norm);
          },
          [w_scale, domain_u_max](double t) -> double {
            long double norm = static_cast<long double>(domain_u_max) / (w_scale + domain_u_max);
            return t / ((w_scale + t) * norm);
          },
          [w_scale, domain_u_max](double l) -> double {
            long double norm = static_cast<long double>(domain_u_max) / (w_scale + domain_u_max);
            return norm * l * w_scale / (1.0 - norm * l);
          },
          domain_u_max};
}

inline warper_product_1d_simple_t make_product_1d_simple_inverse_power(int power, double domain_u_max, double w_scale) {
  if (power == 1) {
    return make_product_1d_simple_inverse(domain_u_max, w_scale);
  }
  auto u_max = static_cast<long double>(domain_u_max);
  return {[power, w_scale, u_max](double t) -> double {
            long double norm =
               (power - 1) * std::pow(u_max + w_scale, power - 1) / (std::pow(u_max / w_scale + 1, power - 1) - 1);
            return norm / std::pow(w_scale + t, power);
          },
          [power, w_scale, u_max](double t) -> double {
            long double norm = std::pow(u_max + w_scale, power - 1) / (std::pow(u_max / w_scale + 1, power - 1) - 1);
            return norm * (std::pow(t / w_scale + 1, power - 1) - 1) / std::pow(t + w_scale, power - 1);
          },
          [power, w_scale, u_max](double l) -> double {
            long double denom =
               std::pow(1 + l * (std::pow(w_scale / (u_max + w_scale), power - 1) - 1), 1.0 / (power - 1.0));
            return w_scale * (1 - denom) / denom;
          },
          domain_u_max};
}

} // namespace keldy::warpers
