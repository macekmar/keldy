/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011-2017 by O. Parcollet
 * Copyright (C) 2018 by Simons Foundation
 *   author : O. Parcollet, P. Dumitrescu, N. Wentzell
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

#include "lapack.h"
#include <complex>

namespace triqs::arrays::lapack::f77 {

  inline void getrs(char trans_a, int N, int NRHS, double *A, int LDA, int *ipiv, double *B, int LDB, int &info) {
    LAPACK_dgetrs(&trans_a, &N, &NRHS, A, &LDA, ipiv, B, &LDB, &info);
  }
  inline void getrs(char trans_a, int N, int NRHS, std::complex<double> *A, int LDA, int *ipiv, std::complex<double> *B, int LDB, int &info) {
    LAPACK_zgetrs(&trans_a, &N, &NRHS, A, &LDA, ipiv, B, &LDB, &info);
  }

} // namespace triqs::arrays::lapack::f77
