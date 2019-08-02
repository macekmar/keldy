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

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <triqs/gfs.hpp>

namespace keldy::details {

class CPP2PY_IGNORE gsl_interp_wrapper_t {

 private:
  gsl_interp_accel *acc;
  gsl_spline *spline;

 public:
  gsl_interp_wrapper_t(const gsl_interp_type *spline_type, const triqs::arrays::array<double, 1> &x_data,
                       const triqs::arrays::array<double, 1> &y_data) {
    // checks on x_data / y_data (contiguous etc.)
    size_t size = x_data.size();
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(spline_type, size);

    gsl_spline_init(spline, x_data.data_start(), y_data.data_start(), size);
  }

  double operator()(double x_in) { return gsl_spline_eval(spline, x_in, acc); }

  // triqs::arrays::array<double, 1> operator()(const triqs::arrays::array<double, 1> &x_in) {
  //   auto F = triqs::arrays::map([&](double x_point) { return gsl_spline_eval(spline, x_point, acc); });
  //   return F(x_in);
  // }

  ~gsl_interp_wrapper_t() {
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
};

} // namespace keldy
