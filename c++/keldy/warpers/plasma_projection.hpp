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

#include "plasma_uv.hpp"
#include "product_1d.hpp"
#include "product_1d_simple.hpp"
#include "../common.hpp"
#include "../interfaces/gsl_interp_wrap.hpp"
#include "../interfaces/gsl_minimize.hpp"
#include "../qrng.hpp"

#include <triqs/gfs.hpp>
#include <triqs/arrays.hpp>
#include <triqs/arrays/mpi.hpp>
#include <mpi/mpi.hpp>
#include <gsl/gsl_min.h>

#include <algorithm>
#include <functional>
#include <numeric>
//#include <iomanip>
#include <complex>
#include <cmath>


namespace keldy::warpers {

using namespace triqs::gfs;
using namespace triqs::arrays;

using gf_t = gf<retime, scalar_real_valued>;

class CPP2PY_IGNORE hist_xi {
 public:
  array<double, 1> bin_times;
  array<double, 1> bin_centers;
  array<double, 1> values;
  array<double, 1> counts; // TODO: should be int, complicates get_xi and others...
  double delta;
  int num_bins;
  int nr_sample_points_warper;

  hist_xi(int num_bins_, int nr_sample_points_warper_)
     : num_bins(num_bins_), nr_sample_points_warper(nr_sample_points_warper_) {
    bin_times = array<double, 1>(num_bins + 1);
    delta = 1.0 / (num_bins - 1);
    for (int i = 0; i < num_bins + 1; i++) {
      bin_times(i) = -delta/2.0 + i*delta;
    }

    bin_centers = array<double, 1>(num_bins);
    for (int i = 0; i < bin_centers.size(); i++) {
      bin_centers(i) = (bin_times(i) + bin_times(i + 1)) / 2.0;
    }

    values = array<double, 1>(num_bins);
    values() = 0;
    counts = array<double, 1>(num_bins); // TODO: should be int
    counts() = 0;
  }
};

// sigma is in unit of time step
auto gaussian = [](double x, double x0, double sigma) { return std::exp(-std::pow((x - x0) / sigma, 2)); };

inline auto CPP2PY_IGNORE gaussian_filter(array<double, 1> &y_in, double sigma) {
  int const N = y_in.size();

  array<double, 1> y_out(N);
  y_out() = 0;
  array<double, 1> norm(N);
  norm() = 0;

  array<double, 1> kernel2(2 * N - 1);
  for (int i = 0; i < 2 * N - 1; ++i) {
    kernel2(i) = gaussian(i, N - 1, sigma);
  }

  norm(0) = sum(kernel2(range(N - 1, 2 * N - 2)));
  for (int i = 1; i < N; ++i) {
    norm(i) = norm(i - 1) + kernel2(N - i - 1) - kernel2(2 * N - i - 1);
  }

  mpi::communicator comm {};
  int mpi_rank = comm.rank();
  int mpi_size = comm.size();

  array<double, 1> kernel(N);
  for (int i = 0; i < N; ++i) {
    if (i % mpi_size == mpi_rank) {
      for (int j = 0; j < N; ++j) {
        kernel(j) = gaussian(j, i, sigma);
      }
      norm(i) = sum(kernel);
      y_out(i) = sum(kernel * y_in);
    }
  }
  y_out = mpi::mpi_all_reduce(y_out, comm);
  norm = mpi::mpi_all_reduce(norm, comm);

  return std::make_pair(y_out, norm);
}

inline array<double, 1> CPP2PY_IGNORE kernel_smoothing(array<double, 1> &y_in, double sigma) {
  auto [y_out, norm] = gaussian_filter(y_in, sigma);
  return y_out / norm;
}

inline double CPP2PY_IGNORE LOOCV(array<double, 1> &y_in, double sigma) {
  auto [y_out, norm] = gaussian_filter(y_in, sigma);

  // We exclude the i-th point
  array<double, 1> y_estim = (y_out - y_in) / (norm - 1.);

  return sum(pow(y_estim - y_in, 2));
}

class warper_plasma_projection_t {
 private:
  int order;
  int nr_sample_points_warper;
  int num_bins;
  double t_max{};
  std::function<std::pair<dcomplex, int>(std::vector<double>)> integrand;
  std::vector<double> fn_integrate_norm{};
  std::vector<gf_t> fn_integrated;
  std::vector<gf_t> fn_integrated_inverse;
  warper_product_1d_simple_t w_warper;
  // warper_train_t warper_train;
  std::vector<hist_xi> xi;
  array<double, 1> values;

  std::vector<gf_t> fn;
  std::vector<double> sigmas;
  mpi::communicator comm {};
 public:

  void gather_data(int npts_mean) {
    int bin;
    double value;
    mpi::communicator comm {};
    int mpi_rank = comm.rank();
    int mpi_size = comm.size();

    sobol g = sobol(order, 0); // TODO: different generators
    for (int i = 0; i < npts_mean; i++) {
      auto w = g();

      if (i % mpi_size != mpi_rank) {
        continue;
      }

      auto u = ui_from_vi(t_max, w_warper.ui_from_li(w));  // TODO: Should have pass warper train
      value = std::abs(integrand(u).first * w_warper.jacobian_reverse(w));
      values(i) = std::real(std::pow(std::complex<double>(0,1), order+1)*integrand(u).first);

      for (int axis = 0; axis < order; axis++) {
        bin = int(floor((w[axis] - (-xi[axis].delta / 2.0)) / xi[axis].delta));
        if (bin < 0 or bin > xi[axis].num_bins - 1) {
          TRIQS_RUNTIME_ERROR << "out of bin range";
        }
        xi[axis].values(bin) += value;
        xi[axis].counts(bin)++;
      }
    }
    for (int axis = 0; axis < order; axis++) {
      xi[axis].values = mpi::mpi_all_reduce(xi[axis].values, comm);
      xi[axis].counts = mpi::mpi_all_reduce(xi[axis].counts, comm);
    }

    for (int axis = 0; axis < order; axis++) {
      for (int bin = 0; bin < xi[axis].num_bins; bin++) {
        if (xi[axis].counts(bin) > 0) {
          xi[axis].values(bin) = xi[axis].values(bin)/xi[axis].counts(bin);
        } else {
          xi[axis].values(bin) = 0.0;
        }
      }
    }
  };

  warper_plasma_projection_t(std::function<std::pair<dcomplex, int>(std::vector<double>)> integrand_,
                             warper_product_1d_simple_t w_warper_, double t_max, int order,
                             int nr_sample_points_warper_, int num_bins, int npts_mean, double sigma,
                             bool optimize_sigma = true)
     : integrand(std::move(integrand_)),
       w_warper(w_warper_),
       t_max(t_max),
       order(order),
       nr_sample_points_warper(nr_sample_points_warper_),
       num_bins(num_bins) {
    // Initialize xi;
    for (int axis = 0; axis < order; axis++) {
      xi.push_back(hist_xi(num_bins, nr_sample_points_warper));
    }
    values = array<double, 1>(npts_mean);
    gather_data(npts_mean);

    // Fill function vectors
    // TODO: should the number of bins be different (smaller) than the number of interpolation points?
    for (int axis = 0; axis < order; axis++) {
      fn_integrated.push_back(gf_t({0.0, 1.0, nr_sample_points_warper}));
      fn_integrate_norm.push_back(0);
      fn_integrated_inverse.push_back(gf_t({0.0, 1.0, 1 + 5 * (nr_sample_points_warper - 1)}));
      fn.push_back(gf_t({0.0, 1.0, nr_sample_points_warper}));
      sigmas.push_back(0);
    }
    update_sigma(sigma, optimize_sigma);

    populate_functions();
  };

  void update_sigma(double sigma, bool optimize_sigma = true) {
    //Smooth and update functions
    for (int axis = 0; axis < order; axis++) {
      sigma /= xi[axis].delta; // change of unit
      if (optimize_sigma) {
        auto f = [&val = xi[axis].values](double s) { return LOOCV(val, s); };

        // TODO: should this be done by all ranks ??
        sigma = details::gsl_minimize(f, 0.5, sigma, 1.0 / xi[axis].delta, 1e-8, 0., 20).x;
        // sigma should not be too small compared to xi[axis].delta, or gaussian is narrower than bin
        if (comm.rank() == 0) {
          std::cout << "Optimal sigma = " << sigma * xi[axis].delta << " (axis " << axis << ")" << std::endl;
        }
      }
      sigmas[axis] = sigma;
    }
  }

  void populate_functions() {

    for (int axis = 0; axis < order; axis++) {
      auto &hist = xi[axis];

      array<double, 1> smoothed_proj;
      smoothed_proj = kernel_smoothing(hist.values, sigmas[axis]);

      // Interpolate convolution
      // linear interp of positive data is always a positive function
      details::gsl_interp_wrapper_t smoothed_proj_interp(gsl_interp_linear, hist.bin_centers, smoothed_proj);

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
      array<double, 1> mesh_time(nr_sample_points_warper);
      for (auto const &[i, t] : itertools::enumerate(fn[axis].mesh())) {
        mesh_time(i) = t;
      }
      details::gsl_interp_wrapper_t interpolate(gsl_interp_akima, data, mesh_time); //gsl_interp_steffen
      for (auto l : fn_integrated_inverse[axis].mesh()) {
        fn_integrated_inverse[axis][l] = interpolate(l);
      }
    }
  };

  std::vector<array<double, 1>> get_xi(int axis) const {
    return {xi[axis].bin_times, xi[axis].values, xi[axis].counts, values};
  };

  std::vector<double> get_sigmas() const {
    std::vector<double> output{};
    for (int i = 0; i < order; ++i) {
      output.emplace_back(sigmas[i] * xi[i].delta);
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
