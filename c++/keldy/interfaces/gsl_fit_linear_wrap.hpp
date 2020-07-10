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

#include <triqs/arrays.hpp>
#include <gsl/gsl_fit.h>
#include <algorithm>

using namespace triqs::arrays;

namespace keldy::details {

struct gsl_fit_linear_result {
  double c0;
  double c1;
  array<double, 2> cov;
  double chisq;
};

/// not used, I kept it for possible future usage
inline gsl_fit_linear_result gsl_fit_wlinear_wrapper(array<double, 1> const &x, array<double, 1> const &y,
                                              array<double, 1> const &w) {
  if (x.size() != y.size() or x.size() != w.size()) {
    TRIQS_RUNTIME_ERROR << "x, y and w should have the same size.";
  }
  size_t n = x.size();
  gsl_fit_linear_result res;
  double cov00, cov01, cov11;

  gsl_fit_wlinear(x.storage().data(), 1, w.storage().data(), 1, y.storage().data(), 1, n, &res.c0, &res.c1, &cov00,
                  &cov01, &cov11, &res.chisq);

  res.cov = array<double, 2>(2, 2);
  res.cov(0, 0) = cov00;
  res.cov(0, 1) = cov01;
  res.cov(1, 0) = cov01;
  res.cov(1, 1) = cov11;

  return res;
};

/// kernel must be odd in size
// in and out must have the same size
inline void local_linear_reg(array<double, 1> const &x, array<double, 1> const &y, array<double, 1> &y_out,
                      array<double, 1> const &kernel) {

  long const N = x.size();
  long const K = kernel.size();
  long const H = K / 2;
  assert(K == 2 * H + 1);
  size_t size_fit;
  size_t start_x;
  size_t start_ker;
  size_t end_ker;
  double c0;
  double c1;
  double cov00, cov01, cov11, chisq;

  for (int i = 0; i < x.size(); ++i) {
    start_ker = std::max(0l, H - i);
    end_ker = std::min(N + H - i - 1, K - 1);
    size_fit = end_ker - start_ker + 1;
    start_x = start_ker - H + i;

    assert(start_x >= 0);
    assert(start_ker >= 0);
    assert(size_fit > 0);
    assert(start_x + size_fit <= N);
    assert(start_ker + size_fit <= K);

    gsl_fit_wlinear(x.storage().data() + start_x, 1, kernel.storage().data() + start_ker, 1,
                    y.storage().data() + start_x, 1, size_fit, &c0, &c1, &cov00, &cov01, &cov11, &chisq);

    y_out(i) = c0 + x(i) * c1;
  }
};

} // namespace keldy::details
