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

#include <gsl/gsl_fit.h>
#include <triqs/arrays.hpp>

using namespace triqs::arrays;

namespace keldy::details {

struct gsl_fit_linear_result {
  double c0;
  double c1;
  array<double, 2> cov;
  double chisq;
};

gsl_fit_linear_result gsl_fit_wlinear_wrapper(array_view<double, 1> const &x, array_view<double, 1> const &y,
                                              array_view<double, 1> const &w) {
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

} // namespace keldy::details
