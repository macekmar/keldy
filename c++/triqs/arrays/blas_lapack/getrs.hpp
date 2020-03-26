/*******************************************************************************
 *
 * TRIQS: A Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (c) 2019 The Simons Foundation
 *   authors: Philipp Dumitrescu
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <complex>
#include <triqs/arrays.hpp>

#include "f77/cxx_interface.hpp"
#include <triqs/arrays/blas_lapack/tools.hpp>
#include <triqs/arrays/blas_lapack/qcache.hpp>

namespace triqs::arrays::lapack {

using namespace blas_lapack_tools;

/**
  * Calls getrs on a matrix or view
  * Takes care of making temporary copies if necessary
  */
template <typename MT>
typename std::enable_if<is_blas_lapack_type<typename MT::value_type>::value, int>::type
getrs(MT &A, MT &B, triqs::arrays::vector<int> &ipiv) {
  // assert(A.memory_layout_is_fortran() && B.memory_layout_is_fortran());  // For now enforce this
  if (A.memory_layout_is_c() || B.memory_layout_is_c())
    TRIQS_RUNTIME_ERROR << "matrices passed to getrs is not in Fortran order";

  reflexive_qcache<MT> Ca(A);
  auto dm = first_dim(Ca());
  if (ipiv.size() != dm)
    TRIQS_RUNTIME_ERROR << "getrs : error in ipiv size : found " << ipiv.size() << " while it should be " << dm;
  reflexive_qcache<MT> Cb(B);
  int info;
  triqs::arrays::lapack::f77::getrs('N', dm, get_n_cols(Cb()), Ca().data_start(), get_ld(Ca()), ipiv.data_start(),
                                    Cb().data_start(), get_ld(Cb()), info);
  return info;
}
} // namespace triqs::arrays::lapack
