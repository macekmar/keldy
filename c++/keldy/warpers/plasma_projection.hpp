/*******************************************************************************
 *
 * keldy
 *
 * Copyright (c) 2020 The Simons Foundation
 * Copyright (c) 2020 CEA: Commissariat à l’énergie atomique
 *                         et aux énergies alternatives
 *   authors: Philipp Dumitrescu, Marjan Macek
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
#include "product_1d.hpp"
#include "../interfaces/gsl_interp_wrap.hpp"
#include <algorithm>
#include <any>
#include <functional>
#include <numeric>
#include <triqs/gfs.hpp>
#include "../qrng.hpp"
#include "plasma_uv.hpp"
#include <iomanip>

namespace keldy::warpers {

using namespace triqs::gfs;

using gf_t = triqs::gfs::gf<triqs::gfs::retime, triqs::gfs::scalar_real_valued>;

struct hist_xi {
  std::vector<double> bins;
  std::vector<double> values;
  std::vector<double> y;
  std::vector<int> counts;
};

inline void bin_values(hist_xi &xi, int axis, std::vector<std::vector<double>> points, std::vector<double> values) {
  bool in_range;
  for (auto [i, point] : itertools::enumerate(points)) {
    in_range = point[axis] >= xi.bins.front() and point[axis] <= xi.bins.back();
    if (in_range) {
      for (auto [j, bin] : itertools::enumerate(xi.values)) {
        if (point[axis] >= xi.bins[j] and point[axis] <= xi.bins[j + 1]) {
          xi.values[j] += values[i];
          xi.counts[j] += 1;
          break;
        }
      }
    } else {
      TRIQS_RUNTIME_ERROR << "out of bin range";
    }
  }
  for (int i = 0; i < xi.values.size(); i++) {
    if (xi.counts[i] != 0) {
      xi.values[i] = xi.values[i] / xi.counts[i];
    } else {
      xi.values[i] = 0;
    }
  }
};

involve void convolve(std::vector<double> &signal, std::vector<double> &result, std::vector<double> &window) {
  std::vector<double> conv;
  double c;
  double norm;
  int Ns = signal.size();
  int Nwl = (window.size() - 1) / 2;
  int Nwr = window.size() - Nwl;
  for (int i = 0; i < Ns; i++) {
    c = 0.0;
    norm = 0.0;
    for (int j = std::max(-Nwl, -Ns + 1 + i); j <= std::min(Nwr, i); j++) {
      c += window[j + Nwl] * signal[i - j];
      norm += window[j + Nwl];
    }
    conv.push_back(c / norm);
  }
  std::copy(conv.begin(), conv.end(), result.begin());
};

class warper_plasma_projection_t {
 private:
  int order;
  int nr_function_sample_points;
  double t_max{};
  std::function<std::pair<dcomplex, int>(std::vector<double>)> integrand;
  std::vector<double> fn_integrate_norm{};
  std::vector<gf_t> fn_integrated;
  std::vector<gf_t> fn_integrated_inverse;
  warper_product_1d_t w_warper;

  std::vector<hist_xi> xi;

  std::vector<gf_t> fn;

 public:
  warper_plasma_projection_t(std::function<std::pair<dcomplex, int>(std::vector<double>)> integrand_, warper_product_1d_t w_warper_,
                             double t_max, int order, int nr_function_sample_points, int npts_mean, int n_window, double sigma)
     : integrand(std::move(integrand_)), w_warper(w_warper_), t_max(t_max), order(order), nr_function_sample_points(nr_function_sample_points) {

    int n_bins = nr_function_sample_points;

    sobol g = sobol(order, 0);

    // Generate points and calculate their values
    std::vector<double> u, w;
    std::vector<std::vector<double>> w_pts;
    std::vector<double> w_vals;
    for (int i = 0; i < npts_mean; i++) {
      w = g();
      w_pts.push_back(w);
      u = ui_from_vi(t_max, w_warper.ui_from_li(w));
      double value = std::abs(integrand(u).first * w_warper.jacobian_reverse(w));
      w_vals.push_back( value );
    }

    // Initialize xi
    for (int axis = 0; axis < order; axis++) {
      fn_integrated.push_back(gf_t({0.0, 1.0, nr_function_sample_points}));

      hist_xi xi_empty;
      double delta = fn_integrated[axis].mesh().delta();
      for (auto const &t : fn_integrated[axis].mesh()) {
        xi_empty.bins.push_back(t - delta / 2.0);
        xi_empty.values.push_back(0);
        xi_empty.y.push_back(0);
        xi_empty.counts.push_back(0);
      }
      xi_empty.bins.push_back(xi_empty.bins.back() + delta);
      xi.push_back(xi_empty);
    }

    // Make histogram
    for (int axis = 0; axis < order; axis++) {
      bin_values(xi[axis], axis, w_pts, w_vals);
    }
    // Do a convolution
    update_sigma(sigma, n_window);
  };

  void update_sigma(double sigma, int n_window) {
    std::cout << std::fixed << std::setw( 15 ) << std::setprecision( 14 );
    // Prepare a Gaussian window for convolution
    std::vector<double> gaussian;
    double x;
    for (int i = 0; i < n_window; i++) {
      x = (-1.0 + 2.0 * i / (n_window - 1)) / sigma;
      gaussian.push_back(std::exp(-std::pow(x, 2.0)));
    }
    

    for (int axis = 0; axis < order; axis++) {
      convolve(xi[axis].values, xi[axis].y, gaussian);

      // -- From here on the code is similar to warper_product_1d_t --

      // Integrate Ansatz using Trapezoid Rule
      fn_integrated.push_back(gf_t({0.0, 1.0, nr_function_sample_points}));
      double delta = fn_integrated[axis].mesh().delta();
      for (auto const &[i, t] : itertools::enumerate(fn_integrated[axis].mesh())) {
        if (i==0) { 
          fn_integrated[axis][t] = 0;
        } else {
          fn_integrated[axis][t] = 0.5*(xi[axis].y[i-1] + xi[axis].y[i])*delta;
        }
      }
      auto &data = fn_integrated[axis].data();

      std::partial_sum(data.begin(), data.end(), data.begin());
      fn_integrate_norm.push_back(data(data.size() - 1));
      data /= fn_integrate_norm[axis]; // normalize each axis separtely

      // Inverse Function via interpolation
      triqs::arrays::array<double, 1> mesh_time(nr_function_sample_points);
      for (auto const &[i, t] : itertools::enumerate(fn_integrated[axis].mesh())) {
        mesh_time(i) = t;
      }

      fn_integrated_inverse.push_back(gf_t({0.0, 1.0, 1 + 5 * (nr_function_sample_points - 1)}));
      details::gsl_interp_wrapper_t interpolate(gsl_interp_akima, data, mesh_time); //gsl_interp_steffen
      for (auto l : fn_integrated_inverse[axis].mesh()) {
        fn_integrated_inverse[axis][l] = interpolate(l);
      }

      // Interpolate projection on CDF
      fn.push_back(gf_t({0.0, 1.0, nr_function_sample_points}));
      for (auto const &[i, t] : itertools::enumerate(fn[axis].mesh())) {
        fn[axis][t] = xi[axis].y[i];
        // std::cout << t << " " << (xi[axis].bins[i] + xi[axis].bins[i+1])/2 << std::endl; // TODO
      }
    }
  };

  std::vector<std::vector<double>> get_xi(int axis) {
    std::vector<double> counts(xi[axis].counts.begin(), xi[axis].counts.end());
    return {xi[axis].bins, xi[axis].values, counts, xi[axis].y};
  };

  std::vector<double> ui_from_li(std::vector<double> const &li_vec) const {
    std::vector<double> result = li_vec;
    for (auto [i, li] : itertools::enumerate(result)) {
      li = fn_integrated_inverse[i](li);
    }
    return result;
  }

  std::vector<double> li_from_ui(std::vector<double> const &ui_vec) const {
    auto result = ui_vec;
    // map ui to vi
    for (auto [i, ui] : itertools::enumerate(result)) {
      ui = fn_integrated[i](ui);
    }
    return result;
  }

  double jacobian_reverse(std::vector<double> const &li_vec) const {
    double result = 1.0;
    for (auto [i, li] : itertools::enumerate(li_vec)) {
      result *= fn_integrate_norm[i] / fn[i](fn_integrated_inverse[i](li));
    }
    return result;
  }

  double operator()(std::vector<double> const &ui_vec) const {
    double result = 1.0;
    // auto vi_vec = vi_from_ui(t_max, ui_vec);
    for (auto [i, vi] : itertools::enumerate(ui_vec)) {
      result *= fn[i](vi);
    }
    return result;
  }
};
} // namespace keldy::warpers
