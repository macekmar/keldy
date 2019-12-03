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

#include "gsl_rng_wrap.hpp"
#include <gsl/gsl_monte_vegas.h>
#include <triqs/utility/first_include.hpp>
#include <vector>

namespace keldy {

struct gsl_monte_vegas_params_wrap {
  double alpha;
  size_t iterations;
  int stage;
  int mode;
  int verbose;
};

class gsl_vegas_wrapper_t {
 private:
  using gsl_integrand_t = double (*)(double *, size_t, void *);

  gsl_integrand_t f_wrap = [](double x[], size_t dim, void *p) -> double {
    auto this_ptr = static_cast<gsl_vegas_wrapper_t *>(p);
    auto v = std::vector(x, x + dim);
    return this_ptr->f_(v);
  };

  std::function<double(std::vector<double>)> f_;

  int dim_;

  gsl_monte_vegas_state *mc_state;
  gsl_rng_wrapper_t rng;
  gsl_monte_function gsl_func;

  std::vector<double> x_lim_lower;
  std::vector<double> x_lim_upper;

  double result = 0.0;
  double error = 0.0;

  uint64_t nr_points_run = 0;

 public:
  gsl_vegas_wrapper_t(std::function<double(std::vector<double>)> f, int dim, double hypercube_extent,
                      std::string const &gsl_rng_name)
     : f_(std::move(f)), dim_(dim), rng(gsl_rng_name), x_lim_lower(dim, 0.0), x_lim_upper(dim, hypercube_extent) {
    mc_state = gsl_monte_vegas_alloc(dim_);

    // Defines Interface for function
    gsl_func.f = f_wrap;
    gsl_func.dim = dim_;
    gsl_func.params = (void *)this;
  }

  double operator()(std::vector<double> x) { return f_(x); }

  double get_result() { return result; }

  double get_error() { return error; }

  uint64_t get_nr_points_run() { return nr_points_run; }

  double chisq() { return gsl_monte_vegas_chisq(mc_state); }

  void run(int nr_steps) {
    gsl_monte_vegas_integrate(&gsl_func, x_lim_lower.data(), x_lim_upper.data(), dim_, nr_steps, rng.get_ptr(),
                              mc_state, &result, &error);
    nr_points_run += nr_steps;
  }

  void set_params(gsl_monte_vegas_params_wrap p) CPP2PY_ARG_AS_DICT {
    auto params_tmp = gsl_monte_vegas_params{p.alpha, p.iterations, p.mode, p.stage, p.verbose, nullptr};
    gsl_monte_vegas_params_set(mc_state, &params_tmp);
  }

  gsl_monte_vegas_params_wrap get_params() {
    gsl_monte_vegas_params p;
    gsl_monte_vegas_params_get(mc_state, &p);
    return gsl_monte_vegas_params_wrap{p.alpha, p.iterations, p.mode, p.stage, p.verbose};
  }

  ~gsl_vegas_wrapper_t() { gsl_monte_vegas_free(mc_state); }
};

} // namespace keldy
