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

#include "common.hpp"
#include "interfaces/gsl_integration_wrap.hpp"
#include <triqs/utility/first_include.hpp>

namespace keldy {
/*
 * Compute the integral along a contour in the complex plane made of three straight lines.
 *
 * The middle segment is [left_turn_pt, right_turn_pt] on the real axis, it is
 * itselfcut at each Fermi level. The two other are semi-infinite lines
 * starting on each end of the middle segment.
 *
 */
class CPP2PY_IGNORE contour_integration_t {

 private:
  details::gsl_integration_wrapper worker;
  double const abstol = 1e-14;
  double const reltol = 1e-8;
  double const left_turn_pt;
  double const right_turn_pt;
  dcomplex result = 0;
  dcomplex abserr = 0;

 public:
  contour_integration_t(double left_turn_pt_, double right_turn_pt_)
     : worker{1000}, left_turn_pt(left_turn_pt_), right_turn_pt(right_turn_pt_) {
    if (not(left_turn_pt < right_turn_pt)) {
      TRIQS_RUNTIME_ERROR << "The order left_turn_pt(" << left_turn_pt << ") < right_turn_pt(" << right_turn_pt
                          << ") should be respected.";
    }
  };

  dcomplex get_result() const { return result; };
  dcomplex get_abserr() const { return abserr; };

  void integrate(std::function<dcomplex(dcomplex)> func, dcomplex left_dir, dcomplex right_dir) {

    // left tail
    auto f1 = [&](double x) -> dcomplex {
      dcomplex omega = left_dir * x + left_turn_pt;
      return -left_dir * func(omega);
    };
    auto [result_1, abserr_1] = worker.qag_si(f1, 0, std::numeric_limits<double>::infinity(), abstol, reltol);

    // middle segment
    auto f2 = [&](double x) -> dcomplex { return func(x); };
    auto [result_2, abserr_2] = worker.qag(f2, left_turn_pt, right_turn_pt, abstol, reltol, GSL_INTEG_GAUSS61);

    // right tail
    auto f3 = [&](double x) -> dcomplex {
      dcomplex omega = right_turn_pt + right_dir * x;
      return right_dir * func(omega);
    };
    auto [result_3, abserr_3] = worker.qag_si(f3, 0, std::numeric_limits<double>::infinity(), abstol, reltol);

    result = result_1 + result_2 + result_3;
    abserr = abserr_1 + abserr_2 + abserr_3;
  }
};

} // namespace keldy
