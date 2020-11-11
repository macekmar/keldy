/*******************************************************************************
 *
 * keldy
 *
 * Copyright (c) 2019-2020 The Simons Foundation
 *   authors: Philipp Dumitrescu
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

#include "integrand_spin_plus_spin_minus_freq.hpp"
#include <triqs/utility/exceptions.hpp>
#include <triqs/arrays.hpp>
#include <triqs/arrays/matrix.hpp>

#include <triqs/arrays/blas_lapack/tools.hpp>
#include <triqs/arrays/blas_lapack/getrs.hpp>

namespace {

inline int GetBit(int in, int offset) { return static_cast<int>((in & (1 << offset)) != 0); }

inline int GetBitParity(unsigned int in) { return 1 - 2 * __builtin_parity(in); }

} // namespace

namespace keldy::impurity_oneband {

/// Returns the first row of the comatrix of `mat`.
// transpose matrix first to get column expansion.
// Strategy: Do row expansion by solving linear equation for cofactors.
// TODO: Do we want to go to full pivoting rather than partial pivoting for LU decomposition?
// TODO: Consider checking Condition Number via LAPACK_zgecon
// TODO: Consider ZGEEQU for equilibiration (scaling of rows / columns for better conditioning)
inline nda::vector<dcomplex> first_row_expansion(nda::matrix<dcomplex> &mat) {
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

  return x_minors(nda::range(), 0);
}

// TODO: Where is sorting of times done?
// TODO: Efficnencies for spin symmetric case.

std::pair<nda::array<dcomplex, 1>, int>
integrand_spin_plus_spin_minus_freq::operator()(std::vector<double> const &times, bool const keep_u_hypercube) const {
  using namespace triqs::arrays;

  nda::array<dcomplex, 1> result(chi.get_nr_omegas());
  result() = 0;

  // Interaction starts a t = 0
  if (keep_u_hypercube) {
    if (std::any_of(times.cbegin(), times.cend(), [](double t) { return t < 0.0; })) {
      return std::make_pair(result, 0);
    }
  }

  int order_n = times.size();

  // TODO: if order_n == 0

  auto a_up = g_idx_X; // This is now a dummy index
  auto a_do = g_idx_X;
  a_up.spin = up;
  a_do.spin = down;
  // define time-splitting for external-points
  a_up.contour.timesplit_n = order_n;
  a_do.contour.timesplit_n = order_n;

  // Pre-Compute Large Matrix.
  matrix<dcomplex> wick_matrix_up(2 * order_n + 1, 2 * order_n + 1);
  matrix<dcomplex> wick_matrix_do(2 * order_n + 1, 2 * order_n + 1);

  // Vector of indices for Green functions
  std::vector<gf_index_t> all_config_up(2 * order_n + 1);
  std::vector<gf_index_t> all_config_do(2 * order_n + 1);
  for (int i = 0; i < order_n; i++) {
    all_config_up[i] = gf_index_t{times[i], up, forward, i};
    all_config_up[i + order_n] = gf_index_t{times[i], up, backward, i};
    all_config_do[i] = gf_index_t{times[i], down, forward, i};
    all_config_do[i + order_n] = gf_index_t{times[i], down, backward, i};
  }
  all_config_up[2 * order_n] = a_up;
  all_config_do[2 * order_n] = a_do;

  // Index for external index in s1
  int external_idx = 2 * order_n;

  // #pragma omp parallel for
  for (int i = 0; i < 2 * order_n + 1; i++) {
    for (int j = 0; j < 2 * order_n + 1; j++) {
      wick_matrix_up(i, j) = g0(all_config_up[i], all_config_up[j]);
      wick_matrix_do(i, j) = g0(all_config_do[i], all_config_do[j]);
    }
  }

  uint64_t nr_keldysh_configs = (uint64_t(1) << order_n);

  std::vector<int> col_pick_up(order_n + 1);
  std::vector<int> col_pick_do(order_n + 1);
  matrix<dcomplex> g_mat_up(order_n + 1, order_n + 1, FORTRAN_LAYOUT);
  matrix<dcomplex> g_mat_do(order_n + 1, order_n + 1, FORTRAN_LAYOUT);

  // Iterate over other Keldysh index configurations. Splict smaller determinant from precomuted matrix
  // #pragma omp parallel for reduce(+: integrand_result)
  for (uint64_t idx_kel = 0; idx_kel < nr_keldysh_configs; idx_kel++) {
    // Indices of Rows / Cols to pick. Cycle through and shift by (0/1) * order_n depending on idx_kel configuration
    col_pick_up[0] = external_idx;
    col_pick_do[0] = external_idx;
    for (int i = 0; i < order_n; i++) {
      col_pick_up[i + 1] = i + GetBit(idx_kel, i) * order_n;
      col_pick_do[i + 1] = i + GetBit(idx_kel, i) * order_n;
    }

    // Extract data into temporary matrices
    for (int i = 0; i < order_n + 1; ++i) {
      for (int j = 0; j < order_n + 1; ++j) {
        g_mat_up(i, j) = wick_matrix_up(col_pick_up[i], col_pick_up[j]);
        g_mat_do(j, i) = wick_matrix_do(col_pick_do[i], col_pick_do[j]); // transposed
      }
    }

    auto cofactors_up = first_row_expansion(g_mat_up);
    auto cofactors_do = first_row_expansion(g_mat_do);

    nda::array<dcomplex, 1> result_tmp(chi.get_nr_omegas());
    result_tmp() = 0;

    int sign = 1;
    for (int k = 0; k < order_n + 1; ++k) {
      sign = -sign;
      auto ind_k = all_config_up[col_pick_up[k]];
      for (int l = 0; l < order_n + 1; ++l) {
        sign = -sign;
        if (k == 0 && l == 0) {
          continue; // vaccum diagrams
        }
        auto ind_l = all_config_do[col_pick_do[l]];
        result_tmp += sign * chi(ind_k, ind_l) * cofactors_up(k) * cofactors_do(l);
      }
    }

    // Multiply by Keldysh index parity
    result_tmp *= GetBitParity(idx_kel);
    result += result_tmp;
  }
  return std::make_pair(result, 1);
}

} // namespace keldy::impurity_oneband
