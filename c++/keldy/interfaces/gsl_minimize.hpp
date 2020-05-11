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

#include <cmath>
#include <keldy/common.hpp>
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>
#include <triqs/utility/exceptions.hpp>
#include <iostream>

namespace keldy::details {

/// Container for the result of gsl_minimize
struct gsl_minimize_result {
  double x = 0;
  double x_lower = 0;
  double x_upper = 0;
  double f = 0;
  double f_lower = 0;
  double f_upper = 0;
  int nr_iter = 0;
  int status = 0;

  friend std::ostream &operator<<(std::ostream &os, gsl_minimize_result const &res) {
    os << "Status: " << res.status << "   Nr iterations: " << res.nr_iter << std::endl;
    os << "Minimum and bounds: x = " << res.x_lower << " < " << res.x << " < " << res.x_upper << std::endl;
    os << "                    f = " << res.f_lower << "   " << res.f << "   " << res.f_upper << std::endl;
    os << std::endl;
    return os;
  };
};

/// Scalar local minimizer using Brent algorithm.
template <typename T,
          typename std::enable_if<std::is_same_v<decltype(std::declval<T>()(std::declval<double>())), double>>::type * =
             nullptr>
gsl_minimize_result gsl_minimize(T const &f, double x_lower, double x_guess, double x_upper, double abstol,
                                 double reltol = 0, int max_iter = 1000, bool end_message = false) {

  gsl_function f_gsl{[](double x, void *param) -> double {
                       auto &f = *static_cast<T *>(param);
                       return f(x);
                     },
                     (void *)&f};

  gsl_min_fminimizer *workspace = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);

  gsl_minimize_result res;
  res.nr_iter = 0;

  auto gsl_previous_handler = gsl_set_error_handler_off(); // temporarily remove error handler

  res.status = gsl_min_fminimizer_set(workspace, &f_gsl, x_guess, x_lower, x_upper);

  if (res.status) {
    if (res.status == GSL_EINVAL) {
      TRIQS_RUNTIME_ERROR << "Endpoints do not enclose a minimum. Make sure f(x_lower) > f(x_guess) < f(x_upper): f("
                          << x_lower << ")=" << f(x_lower) << ", f(" << x_guess << ")=" << f(x_guess) << ", f("
                          << x_upper << ")=" << f(x_upper);
    }
    if (res.status == GSL_EBADFUNC) {
      TRIQS_RUNTIME_ERROR << "Function returned Nan or Inf: f(" << x_lower << ")=" << f(x_lower) << ", f(" << x_guess
                          << ")=" << f(x_guess) << ", f(" << x_upper << ")=" << f(x_upper);
    }
    // unexpected error
    TRIQS_RUNTIME_ERROR << "GSL error " << res.status << " : " << gsl_strerror(res.status);
  };

  do {
    res.nr_iter++;
    res.status = gsl_min_fminimizer_iterate(workspace);

    if (res.status) {
      if (res.status == GSL_EBADFUNC) {
        TRIQS_RUNTIME_ERROR << "Function returned Nan or Inf."; // Don't know where though.
      }
      if (res.status == GSL_FAILURE) {
        TRIQS_RUNTIME_ERROR << "Could not improve the current best approximation or bounding interval.";
      }
      // unexpected error
      TRIQS_RUNTIME_ERROR << "GSL error " << res.status << " : " << gsl_strerror(res.status);
    }

    res.x_lower = gsl_min_fminimizer_x_lower(workspace);
    res.x_upper = gsl_min_fminimizer_x_upper(workspace);

    res.status = gsl_min_test_interval(res.x_lower, res.x_upper, abstol, reltol);

  } while (res.status == GSL_CONTINUE && res.nr_iter < max_iter);

  gsl_set_error_handler(gsl_previous_handler); // put back previous error handler

  if (res.status != GSL_SUCCESS && res.nr_iter != max_iter) {
    TRIQS_RUNTIME_ERROR << "GSL error " << res.status << " : " << gsl_strerror(res.status);
  }

  if (end_message) {
    if (res.status == GSL_SUCCESS) {
      std::cout << "Minimization successful." << std::endl;
    }
    if (res.nr_iter == max_iter) {
      std::cout << "Minimization reached the max number of iterations." << std::endl;
    }
  }

  res.x = gsl_min_fminimizer_x_minimum(workspace);
  res.f = gsl_min_fminimizer_f_minimum(workspace);
  res.f_lower = gsl_min_fminimizer_f_lower(workspace);
  res.f_upper = gsl_min_fminimizer_f_upper(workspace);

  gsl_min_fminimizer_free(workspace);

  return res;
};

} // namespace keldy::details
