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
#include "product_1d.hpp"
#include "../interfaces/gsl_interp_wrap.hpp"
#include <algorithm>
#include <any>
#include <functional>
#include <numeric>
#include <triqs/gfs.hpp>
#include "../qrng.hpp"
#include <triqs/arrays.hpp>

#include <complex>
#include <cmath>

#include <mpi/mpi.hpp>
#include <triqs/arrays/mpi.hpp>

namespace keldy {

using namespace triqs::gfs;

using gf_t = triqs::gfs::gf<triqs::gfs::retime, triqs::gfs::scalar_real_valued>;

class hist_xi {
  public:
  triqs::arrays::array<double, 1> bin_times;
  triqs::arrays::array<double, 1> values;
  triqs::arrays::array<double, 1> counts; // TODO: should be int, complicates get_xi and others...
  triqs::arrays::array<double, 1> y;
  double t_min;
  double t_max;
  double delta;
  int num_bins; 
  hist_xi(double t_min_, double t_max_, int num_bins_) : t_min(t_min_), t_max(t_max_), num_bins(num_bins_) {
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

class warper_plasma_projection_t {
 private:
  int order;
  int num_bins;
  double t_max{};
  std::function<std::pair<dcomplex, int>(std::vector<double>)> integrand;
  std::vector<double> fn_integrate_norm{};
  std::vector<gf_t> fn_integrated;
  std::vector<gf_t> fn_integrated_inverse;
  warper_product_1d_t w_warper;
  // warper_train_t warper_train;
  std::vector<hist_xi> xi;
  triqs::arrays::array<double, 1> values;

  std::vector<gf_t> fn;

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

  warper_plasma_projection_t(std::function<std::pair<dcomplex, int>(std::vector<double>)> integrand_, warper_product_1d_t w_warper_,
                             double t_max, int order, int num_bins, int npts_mean, int n_window_, double sigma)
     : integrand(std::move(integrand_)), w_warper(w_warper_), t_max(t_max), order(order), num_bins(num_bins) {

    int n_window = n_window_ % 2 == 0 ? n_window_ + 1 : n_window_;

    // Initialize xi;
    for (int axis = 0; axis < order; axis++) {
      xi.push_back(hist_xi(0, 1, num_bins));
    }
    values = triqs::arrays::array<double, 1>(npts_mean);
    gather_data(npts_mean);

    // Fill function vectors
    // TODO: should the number of bins be different (smaller) than the number of interpolation points?
    for (int axis = 0; axis < order; axis++) {
      fn_integrated.push_back(gf_t({0.0, 1.0, num_bins}));
      fn_integrate_norm.push_back(0);
      fn_integrated_inverse.push_back(gf_t({0.0, 1.0, 1 + 5 * (num_bins - 1)}));
      fn.push_back(gf_t({0.0, 1.0, num_bins}));
    }
    update_sigma(sigma, n_window);
  };

  void update_sigma(double sigma, int n_window) {   
    // Prepare a window for convolution
    triqs::arrays::array<double, 1> window(n_window);
    auto gaussian_window = [n_window, sigma](int i) { 
      double x = (-1.0 + 2.0 * i / (n_window - 1)) / sigma;
      return std::exp(-std::pow(x, 2.0));};

    for (int i = 0; i < n_window; i++) {
      window(i) = gaussian_window(i);
    }
    
    // Convolve and update functions
    for (int axis = 0; axis < order; axis++) {
      convolve(xi[axis].values, xi[axis].y, window);

      // -- From here on the code is similar to warper_product_1d_t --
      // Integrate Ansatz using Trapezoid Rule
      double delta = fn_integrated[axis].mesh().delta();
      for (auto const &[i, t] : itertools::enumerate(fn_integrated[axis].mesh())) {
        if (i==0) { 
          fn_integrated[axis][t] = 0;
        } else {
          fn_integrated[axis][t] = 0.5*(xi[axis].y(i-1) + xi[axis].y(i))*delta;
        }
      }
      auto &data = fn_integrated[axis].data();

      std::partial_sum(data.begin(), data.end(), data.begin());
      fn_integrate_norm[axis] = data(data.size() - 1);
      data /= fn_integrate_norm[axis]; // normalize each axis separtely

      // Inverse Function via interpolation
      triqs::arrays::array<double, 1> mesh_time(num_bins);
      for (auto const &[i, t] : itertools::enumerate(fn_integrated[axis].mesh())) {
        mesh_time(i) = t;
      }
      
      details::gsl_interp_wrapper_t interpolate(gsl_interp_akima, data, mesh_time); //gsl_interp_steffen
      for (auto l : fn_integrated_inverse[axis].mesh()) {
        fn_integrated_inverse[axis][l] = interpolate(l);
      }

      // Interpolate convolution
      for (auto const &[i, t] : itertools::enumerate(fn[axis].mesh())) {
        fn[axis][t] = xi[axis].y(i);
      }
    }
  };

  std::vector<triqs::arrays::array<double, 1>> get_xi(int axis) {
    return {xi[axis].bin_times, xi[axis].values, xi[axis].counts, xi[axis].y, values};
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

  double jacobian(std::vector<double> const &li_vec) const {
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
} // namespace keldy
