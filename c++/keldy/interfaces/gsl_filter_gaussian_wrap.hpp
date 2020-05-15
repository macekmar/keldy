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

#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_filter.h>
// gsl_filter.h requires GSL 2.6 or higher, previous version are bugged and cannot compile with c++.
// see https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=912341

#include <triqs/arrays.hpp>

int gsl_check_range = 0; // TODO: how does this work?

using namespace triqs::arrays;

namespace keldy::details {

array<double, 1> gsl_vector_2_triqs_array(gsl_vector const &v) {
  array<double, 1> out(v.size);
  for (int i = 0; i < v.size; ++i) {
    out(i) = gsl_vector_get(&v, i);
  }
  return out;
};

gsl_vector *triqs_array_2_gsl_vector(array<double, 1> const &v) {
  gsl_vector *out = gsl_vector_alloc(v.size());
  for (int i = 0; i < v.size(); ++i) {
    gsl_vector_set(out, i, v(i));
  }
  return out;
};

/*
 * size_window is rounded to next odd integer by GSL
 */
array<double, 1> gsl_filter_gaussian_wrapper(array<double, 1> const &x, long size_window, double sigma,
                                             gsl_filter_end_t const endtype, int order = 0) {
  long K = size_window;
  if (K % 2 == 0) {
    K += 1; // this is done by GSL internally anyway, we do it here to get correct alpha
  }
  auto *workspace = gsl_filter_gaussian_alloc(K);
  gsl_vector *x_gsl = triqs_array_2_gsl_vector(x);
  gsl_vector *y_gsl = gsl_vector_alloc(x.size());

  double const alpha = (K - 1) / (2 * sigma);
  gsl_filter_gaussian(endtype, alpha, order, x_gsl, y_gsl, workspace);

  gsl_filter_gaussian_free(workspace);
  gsl_vector_free(x_gsl);

  auto y = gsl_vector_2_triqs_array(*y_gsl);
  gsl_vector_free(y_gsl);

  return y;
};

array<double, 1> gsl_filter_gaussian_kernel_wrapper(long size_window, double sigma, bool normalize, int order = 0) {
  long K = size_window;
  if (K % 2 == 0) {
    K += 1;
  }
  double const alpha = (K - 1) / (2 * sigma);
  gsl_vector *kernel_gsl = gsl_vector_alloc(K);
  gsl_filter_gaussian_kernel(alpha, order, int(normalize), kernel_gsl);

  auto kernel = gsl_vector_2_triqs_array(*kernel_gsl);
  gsl_vector_free(kernel_gsl);
  return kernel;
};

} // namespace keldy::details
