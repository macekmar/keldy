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

#include "warpers_common.hpp"
#include "product_1d_simple.hpp"
#include "../interfaces/gsl_fit_linear_wrap.hpp"
#include "../interfaces/gsl_filter_gaussian_wrap.hpp"
#include "../interfaces/gsl_interp_wrap.hpp"
#include "../interfaces/gsl_minimize.hpp"
#include "../qrng.hpp"
#include "../binner.hpp"

#include <triqs/gfs.hpp>
#include <triqs/arrays.hpp>
#include <triqs/arrays/mpi.hpp>
//#include <mpi/mpi.hpp>

#include <algorithm>
#include <functional>

namespace keldy::warpers {

using namespace triqs::gfs;
using namespace triqs::arrays;

using gf_t = gf<retime, scalar_real_valued>;

auto gaussian = [](double x, double x0, double sigma) { return std::exp(-std::pow((x - x0) / sigma, 2) / 2); };

// sigma is in unit of time step
inline array<double, 1> CPP2PY_IGNORE kernel_smoothing(array<double, 1> const &y_in, double const sigma) {
  int const N = y_in.size();
  array<double, 1> y_out(N);
  array<double, 1> x(N);
  for (int i = 0; i < N; ++i) {
    x(i) = i;
  }

  int const H = std::max(std::min(int(3 * sigma), N), 1);
  int const K = 2 * H + 1;
  array<double, 1> kernel(K);
  for (int j = 0; j < K; ++j) {
    kernel(j) = gaussian(j, H, sigma);
  }

  details::local_linear_reg(x, y_in, y_out, kernel);

  for (int k = 0; k < N; ++k) {
    if (y_out(k) < 0) y_out(k) = 0; // TODO: or should we raise an error?
  }

  return y_out;
}

// sigma is in unit of time step
[[nodiscard]] CPP2PY_IGNORE inline double optimize_sigma(array<double, 1> const &y, double const lower,
                                                         double const sigma, double const upper) {
  int const N = y.size();
  array<double, 1> y_out(N);
  array<double, 1> x(N);
  for (int i = 0; i < N; ++i) {
    x(i) = i;
  }

  std::cout << sigma << std::endl;
  int const H = std::max(std::min(int(3 * sigma), N), 1);
  int const K = 2 * H + 1;
  array<double, 1> kernel(K);

  auto f = [&x, &y, &kernel, &y_out, &H](double sigma) -> double {
    for (int j = 0; j < kernel.size(); ++j) {
      kernel(j) = gaussian(j, H, sigma);
    }
    kernel(H) = 0;
    details::local_linear_reg(x, y, y_out, kernel);

    return sum(pow(y_out - y, 2));
  };

  return details::gsl_minimize(f, lower, sigma, upper, 1e-8, 0., 20).x;
};

class warper_projection_t {
 private:
  int num_bins;
  int order;
  std::vector<double> fn_integrate_norm;
  std::vector<gf_t> fn_integrated;
  std::vector<gf_t> fn_integrated_inverse;
  std::vector<binner::binner_t<1, 0, double>> xi;

  std::vector<gf_t> fn;
  std::vector<double> sigmas;
  mpi::communicator comm{};

 public:
  void CPP2PY_IGNORE gather_data(std::function<dcomplex(std::vector<double>)> const &warped_integrand,
                                 int const nr_samples) {
    double value;
    int mpi_rank = comm.rank();
    int mpi_size = comm.size();

    sobol g = sobol(order, 0); // TODO: different generators
    for (int i = 0; i < nr_samples; i++) {
      auto l = g();

      if (i % mpi_size != mpi_rank) {
        continue;
      }

      value = std::abs(warped_integrand(l));

      for (int axis = 0; axis < order; axis++) {
        xi[axis](l[axis]) << value;
      }
    }
    for (int axis = 0; axis < order; axis++) {
      xi[axis] = mpi_reduce(xi[axis], comm, 0, true);
    }

    for (int axis = 0; axis < order; axis++) {
      auto data = xi[axis].data();             // view
      auto count = xi[axis].nr_values_added(); // view
      for (int bin = 0; bin < xi[axis].get_nr_bins(); bin++) {
        if (count(bin) > 0) {
          data(bin) = data(bin) / count(bin);
        }
      }
    }
  };

  warper_projection_t(std::function<dcomplex(std::vector<double>)> const &warped_integrand, int const order,
                      int const num_bins, int const nr_samples, const double sigma, bool const optimize_sigma = true)
     : num_bins(num_bins), order(order) {
    // Initialize xi;
    double const delta = 1.0 / (num_bins - 1.0);
    for (int axis = 0; axis < order; axis++) {
      xi.push_back(binner::binner_t<1, 0, double>({std::make_tuple(-delta / 2, 1 + delta / 2, num_bins)}));
    }
    gather_data(warped_integrand, nr_samples);

    // Fill function vectors
    for (int axis = 0; axis < order; axis++) {
      fn_integrated.push_back(gf_t({0.0, 1.0, num_bins}));
      fn_integrate_norm.push_back(0);
      fn_integrated_inverse.push_back(gf_t({0.0, 1.0, 1 + 5 * (num_bins - 1)}));
      fn.push_back(gf_t({0.0, 1.0, num_bins}));
    }
    populate_sigmas(sigma, optimize_sigma);

    populate_functions();
  };

  void CPP2PY_IGNORE populate_sigmas(double const sigma, bool const optimize = true) {
    //Smooth and update functions
    for (int axis = 0; axis < order; axis++) {
      sigmas.push_back(sigma / xi[axis].get_bin_size()); // default sigma, change of unit

      if (optimize) {
        // sigma should not be too small compared to xi[axis].delta, or gaussian is narrower than bin
        if (xi[axis].get_bin_size() >= 2) {
          TRIQS_RUNTIME_ERROR << "delta too large";
        }
        sigmas[axis] = optimize_sigma(xi[axis].get_data(), 0.5, sigmas[axis], 1.0 / xi[axis].get_bin_size());

        if (comm.rank() == 0) {
          std::cout << "Optimal sigma = " << sigmas[axis] * xi[axis].get_bin_size() << " (axis " << axis << ")"
                    << std::endl;
        }
      }
    }
  }

  void CPP2PY_IGNORE populate_functions() {

    for (int axis = 0; axis < order; axis++) {

      array<double, 1> smoothed_proj;
      smoothed_proj = kernel_smoothing(xi[axis].get_data(), sigmas[axis]);

      // Interpolate convolution
      // linear interp of positive data is always a positive function
      auto times = xi[axis].get_bin_coord();
      times(0) = 0;                // snap boundaries
      times(times.size() - 1) = 1; // snap boundaries
      details::gsl_interp_wrapper_t smoothed_proj_interp(gsl_interp_linear, times, smoothed_proj);

      for (auto const &t : fn[axis].mesh()) {
        fn[axis][t] = smoothed_proj_interp(t);
      }
      // -- From here on the code is similar to warper_product_1d_t --
      // Integrate Ansatz using Trapezoid Rule
      double delta = fn_integrated[axis].mesh().delta();
      for (auto const &[i, t] : itertools::enumerate(fn_integrated[axis].mesh())) {
        if (i == 0) {
          fn_integrated[axis][t] = 0;
        } else {
          fn_integrated[axis][t] = fn[axis](t - delta / 2) * delta;
        }
      }
      auto &data = fn_integrated[axis].data();

      std::partial_sum(data.begin(), data.end(), data.begin());
      fn_integrate_norm[axis] = data(data.size() - 1);
      data /= fn_integrate_norm[axis]; // normalize each axis separtely

      // Inverse Function via interpolation
      array<double, 1> mesh_time(num_bins);
      for (auto const &[i, t] : itertools::enumerate(fn[axis].mesh())) {
        mesh_time(i) = t;
      }
      details::gsl_interp_wrapper_t interpolate(gsl_interp_akima, data, mesh_time); //gsl_interp_steffen
      for (auto l : fn_integrated_inverse[axis].mesh()) {
        fn_integrated_inverse[axis][l] = interpolate(l);
      }
    }
  };

  auto const &get_xi(int axis) const { return xi[axis]; };

  std::vector<double> get_sigmas() const {
    std::vector<double> output{};
    for (int i = 0; i < order; ++i) {
      // convert sigmas in natural units
      output.emplace_back(sigmas[i] * xi[i].get_bin_size());
    };
    return output;
  }

  auto get_fi(int axis) const { return fn[axis]; };

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

  double jacobian_forward(std::vector<double> const &ui_vec) const {
    double result = 1.0;
    for (auto [i, vi] : itertools::enumerate(ui_vec)) {
      result *= fn[i](vi);
    }
    return result;
  }

  [[nodiscard]] std::pair<std::vector<double>, double> map_reverse(std::vector<double> const &li_vec) const {
    return std::make_pair(ui_from_li(li_vec), jacobian_reverse(li_vec));
  }

  [[nodiscard]] std::pair<std::vector<double>, double> map_forward(std::vector<double> const &ui_vec) const {
    return std::make_pair(li_from_ui(ui_vec), jacobian_forward(ui_vec));
  }
};

} // namespace keldy::warpers
