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

#include <cmath>
#include <keldy/common.hpp>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <triqs/utility/first_include.hpp>
#include <triqs/utility/exceptions.hpp>
#include <limits>
#include <functional>
#include <type_traits>
#include <utility>

namespace keldy::details {

// Propagates GSL error codes into triqs errors
inline void error_handler(int err_code, double result, double abserr) {
  if (err_code != 0) {
    // rounding and maxiter errors should only be warnings
    if (err_code == GSL_EROUND || err_code == GSL_EMAXITER) {
      std::cout << "GSL error " << err_code << " : " << gsl_strerror(err_code) << " (value = " << result
                << ", error = " << abserr << ")" << std::endl;
    } else {
      TRIQS_RUNTIME_ERROR << "GSL error " << err_code << " : " << gsl_strerror(err_code) << " (value = " << result
                          << ", error = " << abserr << ")";
    }
  }
};

class CPP2PY_IGNORE gsl_integration_wrapper {
 private:
  const size_t max_n_intervals_;
  gsl_integration_workspace *wsp;

 public:
  gsl_integration_wrapper(size_t max_n_intervals)
     : max_n_intervals_(max_n_intervals), wsp(gsl_integration_workspace_alloc(max_n_intervals)) {}

  ~gsl_integration_wrapper() { gsl_integration_workspace_free(wsp); }

  // Integrals on finite interval

  template <typename T,
            typename std::enable_if<std::is_same_v<decltype(std::declval<T>()(std::declval<double>())), double>>::type
               * = nullptr>
  [[nodiscard]] std::pair<double, double> qag(T const &f, double lower_limit_a, double upper_limit_b, double epsabs,
                                              double epsrel, int key) {

    gsl_function f_gsl{[](double x, void *param) -> double {
                         auto &f = *static_cast<T *>(param);
                         return f(x);
                       },
                       (void *)&f};

    double result = 0.0;
    double abserr = 0.0;

    auto err_code = gsl_integration_qag(&f_gsl, lower_limit_a, upper_limit_b, epsabs, epsrel, max_n_intervals_, key,
                                        wsp, &result, &abserr);
    error_handler(err_code, result, abserr);
    return std::make_pair(result, abserr);
  }

  template <typename T,
            typename std::enable_if<std::is_same_v<decltype(std::declval<T>()(std::declval<double>())), dcomplex>>::type
               * = nullptr>
  [[nodiscard]] std::pair<dcomplex, dcomplex> qag(T const &f, double lower_limit_a, double upper_limit_b, double epsabs,
                                                  double epsrel, int key) {

    auto f_real = [&f](double t) -> double { return std::real(f(t)); };
    auto [result_real, abserr_real] = qag(f_real, lower_limit_a, upper_limit_b, epsabs, epsrel, key);

    auto f_imag = [&f](double t) -> double { return std::imag(f(t)); };
    auto [result_imag, abserr_imag] = qag(f_imag, lower_limit_a, upper_limit_b, epsabs, epsrel, key);

    return std::make_pair(result_real + 1_j * result_imag, abserr_real + 1_j * abserr_imag);
  }

  // Integrals with singularities or infinite intervals:

  template <typename T,
            typename std::enable_if<std::is_same_v<decltype(std::declval<T>()(std::declval<double>())), double>>::type
               * = nullptr>
  [[nodiscard]] std::pair<double, double> qag_si(T const &f, double lower_limit_a, double upper_limit_b, double epsabs,
                                                 double epsrel) {
    if (lower_limit_a > upper_limit_b) {
      TRIQS_RUNTIME_ERROR << "lower_limit_a should be <= upper_limit_b";
    }

    gsl_function f_gsl{[](double x, void *param) -> double {
                         auto &f = *static_cast<T *>(param);
                         return f(x);
                       },
                       (void *)&f};
    double result = 0.0;
    double abserr = 0.0;
    int err_code;

    bool lower_inf = (lower_limit_a == -std::numeric_limits<double>::infinity());
    bool upper_inf = (upper_limit_b == +std::numeric_limits<double>::infinity());

    // Uses 15-point Gauss-Kronrod for case with inifinate interval (qagi)
    // Uses 21-point Gauss-Kronrod for finite interval [a,b] (qags)
    if (lower_inf) {
      if (upper_inf) {
        err_code = gsl_integration_qagi(&f_gsl, epsabs, epsrel, max_n_intervals_, wsp, &result, &abserr);
      } else {
        err_code =
           gsl_integration_qagil(&f_gsl, upper_limit_b, epsabs, epsrel, max_n_intervals_, wsp, &result, &abserr);
      }
    } else {
      if (upper_inf) {
        err_code =
           gsl_integration_qagiu(&f_gsl, lower_limit_a, epsabs, epsrel, max_n_intervals_, wsp, &result, &abserr);
      } else {
        err_code = gsl_integration_qags(&f_gsl, lower_limit_a, upper_limit_b, epsabs, epsrel, max_n_intervals_, wsp,
                                        &result, &abserr);
      }
    }
    error_handler(err_code, result, abserr);
    return std::make_pair(result, abserr);
  }

  template <typename T,
            typename std::enable_if<std::is_same_v<decltype(std::declval<T>()(std::declval<double>())), dcomplex>>::type
               * = nullptr>
  [[nodiscard]] std::pair<dcomplex, dcomplex> qag_si(T const &f, double lower_limit_a, double upper_limit_b,
                                                     double epsabs, double epsrel) {
    std::function<double(double)> f_real = [&f](double t) { return std::real(f(t)); };
    auto [result_real, abserr_real] = qag_si(f_real, lower_limit_a, upper_limit_b, epsabs, epsrel);

    std::function<double(double)> f_imag = [&f](double t) { return std::imag(f(t)); };
    auto [result_imag, abserr_imag] = qag_si(f_imag, lower_limit_a, upper_limit_b, epsabs, epsrel);

    return std::make_pair(result_real + 1_j * result_imag, abserr_real + 1_j * abserr_imag);
  }
};

} // namespace keldy::gsl