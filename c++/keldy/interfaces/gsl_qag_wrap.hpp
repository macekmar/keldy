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
#include <triqs/utility/first_include.hpp>
#include <functional>

namespace keldy::details {

class CPP2PY_IGNORE gsl_qag_wrapper_t {

 private:
  const size_t limit;
  gsl_integration_workspace *wsp;
  double result = 0;
  double abserr = 0;

 public:
  gsl_qag_wrapper_t(size_t max_intervals_nb)
     : limit(max_intervals_nb), wsp(gsl_integration_workspace_alloc(max_intervals_nb)) {}

  void integrate(std::function<double(double)> const &f, double a, double b, double epsabs, double epsrel, int key) {
    gsl_function f_gsl;
    f_gsl.function = [](double x, void *param) { return (*static_cast<std::function<double(double)> *>(param))(x); };
    f_gsl.params = (void *)&f;

    gsl_integration_qag(&f_gsl, a, b, epsabs, epsrel, limit, key, wsp, &result, &abserr);
  }

  double get_result() { return result; }
  double get_abserr() { return abserr; }

  ~gsl_qag_wrapper_t() { gsl_integration_workspace_free(wsp); }
};

class CPP2PY_IGNORE gsl_qag_cpx_wrapper_t {

 private:
  gsl_qag_wrapper_t worker;
  dcomplex result = 0;
  double abserr_real = 0;
  double abserr_imag = 0;

 public:
  gsl_qag_cpx_wrapper_t(size_t max_intervals_nb) : worker{max_intervals_nb} {};

  void integrate(std::function<dcomplex(double)> const &f, double a, double b, double epsabs, double epsrel, int key) {
    // real part
    auto f_real = [&f](double t) { return std::real(f(t)); };
    worker.integrate(f_real, a, b, epsabs, epsrel, key);
    result = worker.get_result();
    abserr_real = worker.get_abserr();

    // imaginary part
    auto f_imag = [&f](double t) { return std::imag(f(t)); };
    worker.integrate(f_imag, a, b, epsabs, epsrel, key);
    result += 1_j * worker.get_result();
    abserr_imag = worker.get_abserr();
  }

  dcomplex get_result() { return result; }
  double get_abserr_real() { return abserr_real; }
  double get_abserr_imag() { return abserr_imag; }
};

} // namespace keldy::details
