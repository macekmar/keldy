/*******************************************************************************
 *
 * keldy
 *
 * Copyright (c) 2019-2020 The Simons Foundation
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

#include "common.hpp"
#include <triqs/utility/exceptions.hpp>
#include <triqs/arrays.hpp>
#include <triqs/arrays/matrix.hpp>

#include <triqs/arrays/blas_lapack/tools.hpp>
#include <triqs/arrays/blas_lapack/getrs.hpp>

namespace keldy {

namespace nda = triqs::arrays;

/*
 * Returns the first row of the comatrix of `mat`.
 * Get the first column by transposing the matrix first.
 *
 * Strategy: Do row expansion by solving linear equation for cofactors.
 * TODO: Do we want to go to full pivoting rather than partial pivoting for LU decomposition?
 * TODO: Consider checking Condition Number via LAPACK_zgecon
 * TODO: Consider ZGEEQU for equilibiration (scaling of rows / columns for better conditioning)
 */
inline std::pair<nda::vector<dcomplex>, dcomplex> first_row_expansion(nda::matrix<dcomplex> &mat) {
  // Calculate LU decomposition with partial pivoting
  nda::vector<int> pivot_index_array(first_dim(mat));
  int info = nda::lapack::getrf(mat, pivot_index_array, true);
  if (info != 0) {
    TRIQS_RUNTIME_ERROR << "lapack::getrf failed with code " << info;
  }

  // Extract det from LU decompositon (note permutations)
  dcomplex mat_det = 1.0;
  int det_flip = 1;
  for (int i = 0; i < first_dim(mat); i++) {
    mat_det *= mat(i, i);
    if (pivot_index_array(i) != i + 1) { // TODO: check +1 from fortran index
      det_flip = -det_flip;
    }
  }
  mat_det *= det_flip;

  // Find cofactors by solving linear equation Ax = e1 * det(A)
  nda::matrix<dcomplex> x_minors(first_dim(mat), 1, FORTRAN_LAYOUT);
  x_minors() = 0;
  x_minors(0, 0) = mat_det;

  info = nda::lapack::getrs(mat, x_minors, pivot_index_array);
  if (info != 0) {
    TRIQS_RUNTIME_ERROR << "lapack::getrs failed with code " << info;
  }

  return {x_minors(nda::range(), 0), mat_det};
}

} // namespace keldy
