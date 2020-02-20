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
#include <triqs/arrays.hpp>
#include <gsl/gsl_min.h>
#include <complex>
#include <cmath>

#include <mpi/mpi.hpp>
#include <triqs/arrays/mpi.hpp>

namespace keldy::warpers {

using namespace triqs::gfs;

using gf_t = triqs::gfs::gf<triqs::gfs::retime, triqs::gfs::scalar_real_valued>;

class hist_xi {
  public:
  triqs::arrays::array<double, 1> bin_times;
  triqs::arrays::array<double, 1> values;
  triqs::arrays::array<double, 1> counts; // TODO: should be int, complicates get_xi and others...
  triqs::arrays::array<double, 1> y;
  triqs::arrays::array<double, 1> y_interpolated;
  double t_min;
  double t_max;
  double delta;
  int num_bins; 
  int nr_sample_points_warper;
  hist_xi(double t_min_, double t_max_, int num_bins_, int nr_sample_points_warper_) : 
  t_min(t_min_), t_max(t_max_), num_bins(num_bins_), nr_sample_points_warper(nr_sample_points_warper_) {
    bin_times = triqs::arrays::array<double, 1>(num_bins + 1);
    delta = (t_max - t_min)/(num_bins - 1);
    for (int i = 0; i < num_bins + 1; i++) {
      bin_times(i) = -delta/2.0 + i*delta;
    }
    values = triqs::arrays::array<double, 1>(num_bins);
    values() = 0;
    counts = triqs::arrays::array<double, 1>(num_bins); // TODO: should be int
    counts() = 0;
    y = triqs::arrays::array<double, 1>(num_bins);
    y() = 0;
    y_interpolated = triqs::arrays::array<double, 1>(nr_sample_points_warper);
    y_interpolated() = 0;
  }
};


inline void convolve(triqs::arrays::array<double, 1> &signal, triqs::arrays::array<double, 1> &result, triqs::arrays::array<double, 1> &window) {
// TODO: check norm at the right end
  std::vector<double> conv;
  double c;
  double norm;
  int Ns = signal.size();
  int Nwl = (window.size()+1) / 2;
  int Nwr = window.size() - Nwl;
  for (int i = 0; i < Ns; i++) {
    c = 0.0;
    norm = 0.0;
    for (int j = std::max(-Nwl, -Ns + i); j <= std::min(Nwr, i); j++) {
      c += window(j + Nwl) * signal(i - j - 1);
      norm += window(j + Nwl);
    }
    conv.push_back(c / norm);
  }
  std::copy(conv.begin(), conv.end(), result.begin());
};

inline void kernel_smoothing(triqs::arrays::array<double, 1> &x_in, triqs::arrays::array<double, 1> &x_out, triqs::arrays::array<double, 1> &y_in, triqs::arrays::array<double, 1> &y_out, double sigma) {
  auto gaussian = [sigma](double x, double x0) { return std::exp( -std::pow( (x - x0) / sigma, 2 ) ); };
  auto sum = [](double r, double x) {return r + x;};
  triqs::arrays::array<double, 1> kernel(x_in.size());

  for (int i = 0; i < y_out.size(); i++) {
    kernel() = 0.0; 
    for (auto const &[j, x_] : itertools::enumerate(x_in)) {
      kernel(j) = gaussian(x_, x_out(i));
    }
    // \sum_i K(x_0, x_i)*y_i / \sum_i K(x_0, x_i)
    y_out(i) = triqs::arrays::fold(sum)(kernel*y_in, 0) / triqs::arrays::fold(sum)(kernel, 0);
  }
}

inline double LOOCV(triqs::arrays::array<double, 1> x_in, triqs::arrays::array<double, 1> y_in, double sigma) {

  auto gaussian = [sigma](double x, double x0) { return std::exp( -std::pow( (x - x0) / sigma, 2 ) ); };
  auto sum = [](double r, double x) {return r + x;};

  triqs::arrays::array<double, 1> kernel(x_in.size());
  double err = 0;
  for (auto const &[i, y] : itertools::enumerate(y_in)) {
    kernel() = 0.0; 
    for (auto const &[j, x] : itertools::enumerate(x_in)) {
      kernel(j) = gaussian(x, x_in(i));
    }
    kernel(i) = 0; // This excludes i-th point
    auto y_estim = triqs::arrays::fold(sum)(kernel*y_in, 0) / triqs::arrays::fold(sum)(kernel, 0);
    err += std::pow( y - y_estim, 2);
  }

  return err;
}

inline double golden_section(std::function<double(double)> f, double a, double b, int max_iter = 20) {
  //https://en.wikipedia.org/wiki/Golden-section_search
  double gr = 1.618033988749895;
  double c, d;
  double tol = 1e-8;
  int iter = 0;
  c = b - (b - a) / gr;
  d = a + (b - a) / gr;
  while (std::abs(c - d) > tol or iter < max_iter) {
    if (f(c) < f(d)) {
      b = d;
    } else {
      a = c;
    }
    // We recompute both c and d here to avoid loss of precision which may lead to incorrect results or infinite loop
    c = b - (b - a) / gr;
    d = a + (b - a) / gr;
    iter++;
  }

  return (a+b)/2;
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
  triqs::arrays::array<double, 1> values;

  std::vector<gf_t> fn;
  std::vector<double> sigmas;

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
      value = std::abs(integrand(u).first * w_warper.jacobian(w));
      values(i) = std::real(std::pow(std::complex<double>(0,1), order+1)*integrand(u).first);

      for (int axis = 0; axis < order; axis++) {
        bin = int(floor((w[axis] - xi[axis].t_min - (-xi[axis].delta/2.0)) / xi[axis].delta));
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

  warper_plasma_projection_t(std::function<std::pair<dcomplex, int>(std::vector<double>)> integrand_, warper_product_1d_simple_t w_warper_,
                             double t_max, int order, int nr_sample_points_warper_, int num_bins, int npts_mean, double sigma)
     :  integrand(std::move(integrand_)), w_warper(w_warper_), 
        t_max(t_max), order(order), nr_sample_points_warper(nr_sample_points_warper_), 
        num_bins(num_bins) {
    // Initialize xi;
    for (int axis = 0; axis < order; axis++) {
      xi.push_back(hist_xi(0, 1, num_bins, nr_sample_points_warper));
    }
    values = triqs::arrays::array<double, 1>(npts_mean);
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
    update_sigma(sigma, true);
  };

  void update_sigma(double sigma, bool optimize_sigma = true) {   
    mpi::communicator comm {};
    //Smooth and update functions
    for (int axis = 0; axis < order; axis++) {
      triqs::arrays::array<double, 1> bin_centers(xi[axis].y.size());
      bin_centers() = 0;
      for (int i = 0; i < bin_centers.size(); i++) {
        bin_centers(i) = (xi[axis].bin_times(i) + xi[axis].bin_times(i+1))/2.0;
      }

      triqs::arrays::array<double, 1> y(xi[axis].y.size()); // TODO: how to pass xi[axis].values to lambda?
      y() = 0;
      for (int i = 0; i <y.size(); i++) {
       y(i) = xi[axis].values(i);
      }
      auto f = [bin_centers, y](double s){ return LOOCV(bin_centers, y, s);};
      if (optimize_sigma) {
        sigma = golden_section(f, 0.001, 1.0);
      }
      sigmas[axis] = sigma;
      if (comm.rank() == 0) {
      std::cout << "Optimal sigma = " << sigma << std::endl;
      }

      
      triqs::arrays::array<double, 1> mesh_time(nr_sample_points_warper);
      for (auto const &[i, t] : itertools::enumerate(fn[axis].mesh())) {
        mesh_time(i) = t;
      }
      kernel_smoothing(bin_centers, bin_centers, xi[axis].values, xi[axis].y, sigma);
      kernel_smoothing(bin_centers, mesh_time, xi[axis].values, fn[axis].data(), sigma);

      // Interpolate convolution
      // details::gsl_interp_wrapper_t interpolate1(gsl_interp_akima, bin_centers, xi[axis].y);

      // for (auto const &[i, t] : itertools::enumerate(fn[axis].mesh())) {
      //   fn[axis][t] = interpolate1(t);
      // }
      // -- From here on the code is similar to warper_product_1d_t --
      // Integrate Ansatz using Trapezoid Rule
      double delta = fn_integrated[axis].mesh().delta();
      for (auto const &[i, t] : itertools::enumerate(fn_integrated[axis].mesh())) {
        if (i==0) { 
          fn_integrated[axis][t] = 0;
        } else {
          fn_integrated[axis][t] = fn[axis](t - delta / 2) * delta;
        }
      }
      auto &data = fn_integrated[axis].data();
      xi[axis].y_interpolated() = data;

      std::partial_sum(data.begin(), data.end(), data.begin());
      fn_integrate_norm[axis] = data(data.size() - 1);
      data /= fn_integrate_norm[axis]; // normalize each axis separtely

      // Inverse Function via interpolation
      details::gsl_interp_wrapper_t interpolate(gsl_interp_akima, data, mesh_time); //gsl_interp_steffen
      for (auto l : fn_integrated_inverse[axis].mesh()) {
        fn_integrated_inverse[axis][l] = interpolate(l);
      }
    }
  };

  std::vector<triqs::arrays::array<double, 1>> get_xi(int axis) {
    return {xi[axis].bin_times, xi[axis].values, xi[axis].counts, xi[axis].y, xi[axis].y_interpolated, values};
  };
  std::vector<double> get_sigmas() {
    return sigmas;
  }

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
