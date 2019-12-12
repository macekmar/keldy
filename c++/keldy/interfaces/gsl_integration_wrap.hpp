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

#include <keldy/common.hpp>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <triqs/utility/first_include.hpp>
#include <triqs/utility/exceptions.hpp>
#include <functional>

namespace keldy::details {

// Propagates GSL error codes into triqs errors
void error_handler(int err_code, double result, double abserr) {
  if (err_code) {
    // rounding and maxiter errors should only be warnings
    if (err_code == GSL_EROUND || err_code == GSL_EMAXITER) {
      std::cout << "GSL error " << err_code << " : " << gsl_strerror(err_code) << " (value = " << result
                << ", error = " << abserr << ")" << std::endl;
    } else {
      TRIQS_RUNTIME_ERROR << "GSL error " << err_code << " : " << gsl_strerror(err_code) << " (value = " << result
                          << ", error = " << abserr << ")";
    }
  }
}

class CPP2PY_IGNORE gsl_integration_wrapper_t {

 private:
  const size_t max_nb_inter;
  gsl_integration_workspace *wsp;
  double result = 0;
  double abserr = 0;

 public:
  gsl_integration_wrapper_t(size_t max_nb_inter_)
     : max_nb_inter(max_nb_inter_), wsp(gsl_integration_workspace_alloc(max_nb_inter_)) {}
  ~gsl_integration_wrapper_t() { gsl_integration_workspace_free(wsp); }

  double get_result() { return result; }
  double get_abserr() { return abserr; }

  void integrate_qag(std::function<double(double)> const &f, double a, double b, double epsabs, double epsrel,
                     int key) {
    gsl_function f_gsl;
    f_gsl.function = [](double x, void *param) { return (*static_cast<std::function<double(double)> *>(param))(x); };
    f_gsl.params = (void *)&f;
    auto err_code = gsl_integration_qag(&f_gsl, a, b, epsabs, epsrel, max_nb_inter, key, wsp, &result, &abserr);
    error_handler(err_code, result, abserr);
  }

  void integrate_qagi(std::function<double(double)> const &f, double epsabs, double epsrel) {
    gsl_function f_gsl;
    f_gsl.function = [](double x, void *param) { return (*static_cast<std::function<double(double)> *>(param))(x); };
    f_gsl.params = (void *)&f;
    auto err_code = gsl_integration_qagi(&f_gsl, epsabs, epsrel, max_nb_inter, wsp, &result, &abserr);
    error_handler(err_code, result, abserr);
  }

  void integrate_qagiu(std::function<double(double)> const &f, double a, double epsabs, double epsrel) {
    gsl_function f_gsl;
    f_gsl.function = [](double x, void *param) { return (*static_cast<std::function<double(double)> *>(param))(x); };
    f_gsl.params = (void *)&f;
    auto err_code = gsl_integration_qagiu(&f_gsl, a, epsabs, epsrel, max_nb_inter, wsp, &result, &abserr);
    error_handler(err_code, result, abserr);
  }
};

// Integration of complex valued functions
class CPP2PY_IGNORE gsl_integration_cpx_wrapper_t {

 private:
  gsl_integration_wrapper_t worker;
  dcomplex result = 0;
  double abserr_real = 0;
  double abserr_imag = 0;

 public:
  gsl_integration_cpx_wrapper_t(size_t max_intervals_nb) : worker{max_intervals_nb} {};

  dcomplex get_result() { return result; }
  double get_abserr_real() { return abserr_real; }
  double get_abserr_imag() { return abserr_imag; }

  void integrate_qag(std::function<dcomplex(double)> const &f, double a, double b, double epsabs, double epsrel,
                     int key) {
    // real part
    auto f_real = [&f](double t) { return std::real(f(t)); };
    worker.integrate_qag(f_real, a, b, epsabs, epsrel, key);
    result = worker.get_result();
    abserr_real = worker.get_abserr();

    // imaginary part
    auto f_imag = [&f](double t) { return std::imag(f(t)); };
    worker.integrate_qag(f_imag, a, b, epsabs, epsrel, key);
    result += 1_j * worker.get_result();
    abserr_imag = worker.get_abserr();
  }

  void integrate_qagi(std::function<dcomplex(double)> const &f, double epsabs, double epsrel) {
    // real part
    auto f_real = [&f](double t) { return std::real(f(t)); };
    worker.integrate_qagi(f_real, epsabs, epsrel);
    result = worker.get_result();
    abserr_real = worker.get_abserr();

    // imaginary part
    auto f_imag = [&f](double t) { return std::imag(f(t)); };
    worker.integrate_qagi(f_imag, epsabs, epsrel);
    result += 1_j * worker.get_result();
    abserr_imag = worker.get_abserr();
  }

  void integrate_qagiu(std::function<dcomplex(double)> const &f, double a, double epsabs, double epsrel) {
    // real part
    auto f_real = [&f](double t) { return std::real(f(t)); };
    worker.integrate_qagiu(f_real, a, epsabs, epsrel);
    result = worker.get_result();
    abserr_real = worker.get_abserr();

    // imaginary part
    auto f_imag = [&f](double t) { return std::imag(f(t)); };
    worker.integrate_qagiu(f_imag, a, epsabs, epsrel);
    result += 1_j * worker.get_result();
    abserr_imag = worker.get_abserr();
  }
};

} // namespace keldy::details
