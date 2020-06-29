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
#include <triqs/utility/exceptions.hpp>
#include <triqs/utility/macros.hpp>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

namespace keldy::warpers {

using gf_t = triqs::gfs::gf<triqs::gfs::retime, triqs::gfs::scalar_real_valued>;

class warper_product_1d_simple_t {
 protected:
  double f1_integrated_normalization = 1.0;

  std::function<double(double)> f1 = [](double /*u*/) { return 1.0; };
  std::function<double(double)> f1_integrated = [](double u) { return u; };
  std::function<double(double)> f1_integrated_inverse = [](double l) { return l; };

  constexpr static double domain_u_min = 0.0; // alwyas start at 0.0
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
      xi = f1_integrated_inverse(xi);
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

class warper_product_1d_simple_interp_nearest_t : public warper_product_1d_simple_t {
 private:
  int nr_sample_points = 0;

  nda::vector<double> times_u_pts{};
  nda::vector<double> f1_pts{};
  nda::vector<double> f1_integrated_pts{};

  // f is constant interpolated to nearest downward: this piggy-backs on the acc lookup for f_integrated

  std::unique_ptr<gsl_spline, decltype(&gsl_spline_free)> sp_f_integrated;
  std::unique_ptr<gsl_spline, decltype(&gsl_spline_free)> sp_f_integrated_inverse;

  std::unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)> acc_f_integrated;
  std::unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)> acc_f_integrated_inverse;

 public:
  warper_product_1d_simple_interp_nearest_t()
     : sp_f_integrated{gsl_spline_alloc(gsl_interp_linear, nr_sample_points), gsl_spline_free},
       sp_f_integrated_inverse{gsl_spline_alloc(gsl_interp_linear, nr_sample_points), gsl_spline_free},
       acc_f_integrated{gsl_interp_accel_alloc(), gsl_interp_accel_free},
       acc_f_integrated_inverse{gsl_interp_accel_alloc(), gsl_interp_accel_free} {};

  warper_product_1d_simple_interp_nearest_t(std::function<double(double)> const &f1_, double domain_u_max_,
                                            int nr_sample_points_)
     : warper_product_1d_simple_t{
        [this](double u) {
          auto i = gsl_interp_accel_find(acc_f_integrated.get(), times_u_pts.data_start(), nr_sample_points, u);
          if (i == times_u_pts.size() - 1) {
            return f1_pts[i];
          };
          return 0.5 * (f1_pts[i] + f1_pts[i + 1]);
        },
        [this](double u) { return gsl_spline_eval(sp_f_integrated.get(), u, acc_f_integrated.get()); },
        [this](double l) { return gsl_spline_eval(sp_f_integrated_inverse.get(), l, acc_f_integrated_inverse.get()); },
        domain_u_max_, false},
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
      f1_pts[i] = f1_(times_u_pts[i]);
    }

    f1_integrated_pts[0] = 0.0; // scale co-ordinates later
    for (int i = 1; i < nr_sample_points; i++) {
      f1_integrated_pts[i] =
         f1_integrated_pts[i - 1] + 0.5 * (f1_pts[i - 1] + f1_pts[i]) * (times_u_pts[i] - times_u_pts[i - 1]);
    }

    f1_integrated_normalization = f1_integrated_pts[nr_sample_points - 1];

    //#pragma omp parallel for
    for (int i = 0; i < nr_sample_points; i++) {
      f1_integrated_pts[i] *= codomain_l_max / f1_integrated_normalization;
      f1_pts[i] *= codomain_l_max / f1_integrated_normalization;
    }

    gsl_spline_init(sp_f_integrated.get(), times_u_pts.data_start(), f1_integrated_pts.data_start(), nr_sample_points);
    gsl_spline_init(sp_f_integrated_inverse.get(), f1_integrated_pts.data_start(), times_u_pts.data_start(),
                    nr_sample_points);
  }

  // Need to deep-copy for python layer
  warper_product_1d_simple_interp_nearest_t(const warper_product_1d_simple_interp_nearest_t &o)
     : warper_product_1d_simple_t{
        [this](double u) {
          auto i = gsl_interp_accel_find(acc_f_integrated.get(), times_u_pts.data_start(), nr_sample_points, u);
          if (i == times_u_pts.size() - 1) {
            return f1_pts[i];
          };
          return 0.5 * (f1_pts[i] + f1_pts[i + 1]);
        },
        [this](double u) { return gsl_spline_eval(sp_f_integrated.get(), u, acc_f_integrated.get()); },
        [this](double l) { return gsl_spline_eval(sp_f_integrated_inverse.get(), l, acc_f_integrated_inverse.get()); },
        o.domain_u_max, false},
       nr_sample_points{o.nr_sample_points},
       times_u_pts(o.times_u_pts),
       f1_pts(o.f1_pts),
       f1_integrated_pts(o.f1_integrated_pts),
       sp_f_integrated{gsl_spline_alloc(gsl_interp_linear, nr_sample_points), gsl_spline_free},
       sp_f_integrated_inverse{gsl_spline_alloc(gsl_interp_linear, nr_sample_points), gsl_spline_free},
       acc_f_integrated{gsl_interp_accel_alloc(), gsl_interp_accel_free},
       acc_f_integrated_inverse{gsl_interp_accel_alloc(), gsl_interp_accel_free} {
    gsl_spline_init(sp_f_integrated.get(), times_u_pts.data_start(), f1_integrated_pts.data_start(), nr_sample_points);
    gsl_spline_init(sp_f_integrated_inverse.get(), f1_integrated_pts.data_start(), times_u_pts.data_start(),
                    nr_sample_points);
  }

  warper_product_1d_simple_interp_nearest_t &operator=(const warper_product_1d_simple_interp_nearest_t &o) {
    warper_product_1d_simple_interp_nearest_t tmp(o);
    std::swap(*this, tmp);
    return *this;
  }
};

class warper_product_1d_simple_interp_hybrid_t : public warper_product_1d_simple_t {
 private:
  int nr_sample_points = 0;

  nda::vector<double> times_u_pts{};
  nda::vector<double> f1_integrated_pts{};

  std::unique_ptr<gsl_spline, decltype(&gsl_spline_free)> sp_f_integrated;
  std::unique_ptr<gsl_spline, decltype(&gsl_spline_free)> sp_f_integrated_inverse;

  std::unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)> acc_f_integrated;
  std::unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)> acc_f_integrated_inverse;

 public:
  warper_product_1d_simple_interp_hybrid_t()
     : sp_f_integrated{gsl_spline_alloc(gsl_interp_steffen, nr_sample_points), gsl_spline_free},
       sp_f_integrated_inverse{gsl_spline_alloc(gsl_interp_steffen, nr_sample_points), gsl_spline_free},
       acc_f_integrated{gsl_interp_accel_alloc(), gsl_interp_accel_free},
       acc_f_integrated_inverse{gsl_interp_accel_alloc(), gsl_interp_accel_free} {};

  warper_product_1d_simple_interp_hybrid_t(std::function<double(double)> const &f1_, double domain_u_max_,
                                           int nr_sample_points_)
     : warper_product_1d_simple_t{
        {}, // set f1 later due to normalization
        [this](double u) { return gsl_spline_eval(sp_f_integrated.get(), u, acc_f_integrated.get()); },
        [this](double l) { return gsl_spline_eval(sp_f_integrated_inverse.get(), l, acc_f_integrated_inverse.get()) },
        domain_u_max_,
        false},
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

    // Inetegrate function using Simpsons rule: Evaluate f1 two times for each output point
    f1_integrated_pts[0] = 0.0;      // scale co-ordinates later
    double f2 = f1_(times_u_pts[0]); // temp -- carry over
    for (int i = 1; i < nr_sample_points; i++) {
      double f0 = f2;
      double f1 = f1_(0.5 * (times_u_pts[i - 1] + times_u_pts[i]));
      f2 = f1_(times_u_pts[i]);
      f1_integrated_pts[i] = f1_integrated_pts[i - 1] + delta / 6. * (f0 + 4 * f1 + f2);
    }

    f1_integrated_normalization = f1_integrated_pts[nr_sample_points - 1];

    //#pragma omp parallel for
    for (int i = 0; i < nr_sample_points; i++) {
      f1_integrated_pts[i] = (f1_integrated_pts[i] / f1_integrated_normalization) * codomain_l_max;
    }

    // Ensure that waf1_integrated_pts[i] == 1.0

    long double total_normalization = codomain_l_max / f1_integrated_normalization;
    f1 = [total_normalization, f1_](double u) { return f1_(u) * total_normalization; };

    gsl_spline_init(sp_f_integrated.get(), times_u_pts.data_start(), f1_integrated_pts.data_start(), nr_sample_points);
    gsl_spline_init(sp_f_integrated_inverse.get(), f1_integrated_pts.data_start(), times_u_pts.data_start(),
                    nr_sample_points);
  }

  // Need to deep-copy for python layer
  warper_product_1d_simple_interp_hybrid_t(const warper_product_1d_simple_interp_hybrid_t &o)
     : warper_product_1d_simple_t{
        o.f1, [this](double u) { return gsl_spline_eval(sp_f_integrated.get(), u, acc_f_integrated.get()); },
        [this](double l) { return gsl_spline_eval(sp_f_integrated_inverse.get(), l, acc_f_integrated_inverse.get()); },
        o.domain_u_max, false},
       nr_sample_points{o.nr_sample_points},
       times_u_pts(o.times_u_pts),
       f1_integrated_pts(o.f1_integrated_pts),
       sp_f_integrated{gsl_spline_alloc(gsl_interp_steffen, nr_sample_points), gsl_spline_free},
       sp_f_integrated_inverse{gsl_spline_alloc(gsl_interp_steffen, nr_sample_points), gsl_spline_free},
       acc_f_integrated{gsl_interp_accel_alloc(), gsl_interp_accel_free},
       acc_f_integrated_inverse{gsl_interp_accel_alloc(), gsl_interp_accel_free} {
    gsl_spline_init(sp_f_integrated.get(), times_u_pts.data_start(), f1_integrated_pts.data_start(), nr_sample_points);
    gsl_spline_init(sp_f_integrated_inverse.get(), f1_integrated_pts.data_start(), times_u_pts.data_start(),
                    nr_sample_points);
  }

  warper_product_1d_simple_interp_hybrid_t &operator=(const warper_product_1d_simple_interp_hybrid_t &o) {
    warper_product_1d_simple_interp_hybrid_t tmp(o);
    std::swap(*this, tmp);
    return *this;
  }
};

// ***************************************************************************************************************

// Maker Functions:

inline warper_product_1d_simple_t make_product_1d_simple_exponential_nointerp(double domain_u_max, double w_scale) {
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

inline warper_product_1d_simple_t make_product_1d_simple_inverse_nointerp(double domain_u_max, double w_scale) {
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

inline warper_product_1d_simple_t make_product_1d_simple_inverse_square_nointerp(double domain_u_max, double w_scale) {
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

} // namespace keldy::warpers
