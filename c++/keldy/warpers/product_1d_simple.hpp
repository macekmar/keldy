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
#include "../interfaces/gsl_interp_wrap.hpp"
#include <algorithm>
#include <any>
#include <functional>
#include <numeric>
#include <triqs/gfs.hpp>

namespace keldy::warpers {

using gf_t = triqs::gfs::gf<triqs::gfs::retime, triqs::gfs::scalar_real_valued>;

class warper_product_1d_simple_t {
 private:
  double t_max{};
  double f1_integrated_norm{};
  gf_t f1_integrated;
  gf_t f1_integrated_inverse;
  std::function<double(double)> f1;

 public:
  // Identity Constructor: should use this if nothing else is specified
  warper_product_1d_simple_t(double t_max_)
     : warper_product_1d_simple_t{identity_function{}, linear_function{}, linear_function{}, t_max_, 4} {}

  // Pass Warping Function: numerically perform integration
  // TODO: points vs resampling points
  warper_product_1d_simple_t(std::function<double(double)> f1_, double t_max_, int nr_function_sample_points)
     : t_max(t_max_), f1(std::move(f1_)) {

    // Integrate Ansatz using Trapezoid Rule
    f1_integrated = gf_t({0.0, t_max, nr_function_sample_points});
    double delta = f1_integrated.mesh().delta();
    for (auto const &t : f1_integrated.mesh()) {
      f1_integrated[t] = f1(t - delta / 2) * delta;
    }

    auto &data = f1_integrated.data();
    data(0) = 0.0;

    std::partial_sum(data.begin(), data.end(), data.begin());
    f1_integrated_norm = data(data.size() - 1);
    data /= f1_integrated_norm; // normalize

    // Inverse Function via interpolation
    triqs::arrays::array<double, 1> mesh_time(nr_function_sample_points);
    for (auto const &[i, t] : itertools::enumerate(f1_integrated.mesh())) {
      mesh_time(i) = t;
    }

    f1_integrated_inverse = gf_t({0.0, 1.0, 1 + 5 * (nr_function_sample_points - 1)});
    details::gsl_interp_wrapper_t interpolate(gsl_interp_akima, data, mesh_time); //gsl_interp_steffen
    for (auto l : f1_integrated_inverse.mesh()) {
      f1_integrated_inverse[l] = interpolate(l);
    }
  }

  // Constructor to use if f1, f1_integrated and f1_integrated_inverse can be provided analytically
  // To do:points vs resampling points
  warper_product_1d_simple_t(std::function<double(double)> f1_, std::function<double(double)> f1_integrated_,
                             std::function<double(double)> f1_integrated_inverse_, double t_max_,
                             int nr_function_sample_points)
     : t_max(t_max_), f1_integrated_norm(f1_integrated_(t_max)), f1(std::move(f1_)) {

    // check consistency of functions provided
    double l_prime = 0.;
    for (double l : {0., 0.25, 0.5, 0.75, 1.}) {
      l_prime = f1_integrated_norm * l + f1_integrated_(0) * (1 - l);
      if (std::abs(f1_integrated_(f1_integrated_inverse_(l_prime)) - l_prime) > 1e-10) {
        TRIQS_RUNTIME_ERROR << "Inconsistent functions: f1_integrated_inverse should be the inverse of f1_integrated ("
                            << f1_integrated_(f1_integrated_inverse_(l_prime)) << " != " << l_prime << ")";
      }
    }

    f1_integrated = gf_t({0.0, t_max, nr_function_sample_points});
    for (auto const &t : f1_integrated.mesh()) {
      f1_integrated[t] = (f1_integrated_(t) - f1_integrated_(0.)) / (f1_integrated_norm - f1_integrated_(0.));
    }

    f1_integrated_inverse = gf_t({0.0, 1.0, 1 + 5 * (nr_function_sample_points - 1)});
    for (auto const &l : f1_integrated_inverse.mesh()) {
      f1_integrated_inverse[l] = f1_integrated_inverse_(f1_integrated_norm * l + f1_integrated(0) * (1 - l));
    }
  }

  [[nodiscard]] std::pair<std::vector<double>, double> map_reverse(std::vector<double> const &li_vec) const {
    return std::make_pair(ui_from_li(li_vec), jacobian_reverse(li_vec));
  }

  [[nodiscard]] std::pair<std::vector<double>, double> map_forward(std::vector<double> const &ui_vec) const {
    return std::make_pair(li_from_ui(ui_vec), jacobian_forward(ui_vec));
  }

  std::vector<double> ui_from_li(std::vector<double> const &li_vec) const {
    std::vector<double> result = li_vec;
    for (auto &li : result) {
      li = f1_integrated_inverse(li);
    }
    return result;
  }

  std::vector<double> li_from_ui(std::vector<double> const &ui_vec) const {
    auto result = ui_vec;
    // map ui to vi
    for (auto &ui : result) {
      ui = f1_integrated(ui);
    }
    return result;
  }

  double jacobian_reverse(std::vector<double> const &li_vec) const {
    double result = 1.0;
    for (auto li : li_vec) {
      result *= f1_integrated_norm / f1(f1_integrated_inverse(li));
    }
    return result;
  }

  double jacobian_forward(std::vector<double> const &ui_vec) const {
    double result = 1.0;
    for (auto &vi : ui_vec) {
      result *= f1(vi);
    }
    return result;
  }
};

// Maker Functions:

inline warper_product_1d_simple_t make_product_1d_simple_exponential(double time, double w_scale,
                                                                     int nr_sample_points_warper) {
  return {[w_scale](double t) -> double { return std::exp(-(t / w_scale)); },
          [w_scale](double t) -> double { return w_scale * (1 - std::exp(-t / w_scale)); },
          [w_scale](double l) -> double { return -w_scale * std::log(1 - l / w_scale); }, time,
          nr_sample_points_warper};
}

inline warper_product_1d_simple_t make_product_1d_simple_inverse(double time, double w_scale,
                                                                 int nr_sample_points_warper) {
  return {[w_scale](double t) -> double { return w_scale / (w_scale + t); },
          [w_scale](double t) -> double { return w_scale * std::log(1. + t / w_scale); },
          [w_scale](double l) -> double { return w_scale * (std::exp(l / w_scale) - 1.); }, time,
          nr_sample_points_warper};
}

inline warper_product_1d_simple_t make_product_1d_simple_inverse_square(double time, double w_scale,
                                                                        int nr_sample_points_warper) {
  return {[w_scale](double t) -> double { return w_scale * w_scale / ((w_scale + t) * (w_scale + t)); },
          [w_scale](double t) -> double { return w_scale * t / (w_scale + t); },
          [w_scale](double l) -> double { return w_scale * l / (w_scale - l); }, time, nr_sample_points_warper};
}

// (std::string const &label, integrand_g_direct const &f, double time, int nr_sample_points_warper, double w_scale) {
//   if (label == "first_order") {
//     return {[time, &f](double t) -> double { return std::abs(f(std::vector<double>{time - t}).first) + 1e-12; }, time,
//             nr_sample_points_warper};
//   }
// }

} // namespace keldy::warpers
