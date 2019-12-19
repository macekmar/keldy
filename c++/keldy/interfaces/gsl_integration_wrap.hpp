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
#include <functional>

namespace keldy::details {

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

  void integrate_qag(std::function<double(double)> const &f, double a, double b, double epsabs, double epsrel, int key);
  void integrate_qagi(std::function<double(double)> const &f, double epsabs, double epsrel);
  void integrate_qagiu(std::function<double(double)> const &f, double a, double epsabs, double epsrel);
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
                     int key);
  void integrate_qagi(std::function<dcomplex(double)> const &f, double epsabs, double epsrel);
  void integrate_qagiu(std::function<dcomplex(double)> const &f, double a, double epsabs, double epsrel);
};

} // namespace keldy::details
