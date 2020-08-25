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

#include "wick_kernel_single_omega.hpp"
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

// TODO: Where is sorting of times done?
// TODO: Efficnencies for spin symmetric case.

std::pair<dcomplex, int> integrand_g_kernel_single_omega::operator()(std::vector<double> const &times,
                                                                     bool const keep_u_hypercube) const {
  using namespace triqs::arrays;

  // Interaction starts a t = 0
  if (keep_u_hypercube) {
    if (std::any_of(times.cbegin(), times.cend(), [](double t) { return t < 0.0; })) {
      return std::make_pair(dcomplex{0}, 0);
    }
  }

  int order_n = times.size();

  // TODO: if order_n == 0

  dcomplex result = 0;

  auto a = g_idx_X; // This is now a dummy index
  auto b = g_idx_X;
  // define time-splitting for external-points
  a.contour.timesplit_n = order_n;
  b.contour.timesplit_n = order_n;

  // Pre-Comute Large Matrix.
  // "s1": Same spin as external indices / "s2": Opposite spin
  matrix<dcomplex> wick_matrix_s1(2 * order_n + 1, 2 * order_n + 1);
  matrix<dcomplex> wick_matrix_s2(2 * order_n, 2 * order_n);

  // Vector of indices for Green functions
  std::vector<gf_index_t> all_config_1(2 * order_n);
  std::vector<gf_index_t> all_config_2(2 * order_n);
  for (int i = 0; i < order_n; i++) {
    all_config_1[i] = gf_index_t{times[i], a.spin, forward, i};
    all_config_1[i + order_n] = gf_index_t{times[i], a.spin, backward, i};
    all_config_2[i] = gf_index_t{times[i], spin_t(1 - a.spin), forward, i};
    all_config_2[i + order_n] = gf_index_t{times[i], spin_t(1 - a.spin), backward, i};
  }

  // Index for external index in s1
  int external_idx = 2 * order_n;

  wick_matrix_s1(external_idx, external_idx) = g0(a, b, false);
  // #pragma omp parallel for
  for (int i = 0; i < 2 * order_n; i++) {
    wick_matrix_s1(external_idx, i) = g0(a, all_config_1[i]);
    wick_matrix_s1(i, external_idx) = g0(all_config_1[i], b);
    for (int j = 0; j < 2 * order_n; j++) {
      wick_matrix_s1(i, j) = g0(all_config_1[i], all_config_1[j]);
      wick_matrix_s2(i, j) = g0(all_config_2[i], all_config_2[j]);
    }
  }

  uint64_t nr_keldysh_configs = (uint64_t(1) << order_n);

  // Iterate over other Keldysh index configurations. Splict smaller determinant from precomuted matrix
  // #pragma omp parallel for reduce(+: integrand_result)
  for (uint64_t idx_kel = 0; idx_kel < nr_keldysh_configs; idx_kel++) {
    // Indices of Rows / Cols to pick. Cycle through and shift by (0/1) * order_n depending on idx_kel configuration
    std::vector<int> col_pick_s1(order_n + 1);
    std::vector<int> col_pick_s2(order_n);

    col_pick_s1[0] = external_idx;
    for (int i = 0; i < order_n; i++) {
      col_pick_s1[i + 1] = i + GetBit(idx_kel, i) * order_n;
      col_pick_s2[i] = col_pick_s1[i + 1];
    }

    // Extract data into temporary matrices
    matrix<dcomplex> g_mat_s1(order_n + 1, order_n + 1, FORTRAN_LAYOUT);
    matrix<dcomplex> g_mat_s2(order_n, order_n, FORTRAN_LAYOUT);

    for (int i = 0; i < order_n; ++i) {
      for (int j = 0; j < order_n; ++j) {
        g_mat_s2(i, j) = wick_matrix_s2(col_pick_s2[i], col_pick_s2[j]);
      }
    }

    for (int i = 0; i < order_n + 1; ++i) {
      for (int j = 0; j < order_n + 1; ++j) {
        g_mat_s1(i, j) = wick_matrix_s1(col_pick_s1[i], col_pick_s1[j]);
      }
    }

    // Strategy: Do row expansion by solving linear equation for cofactors.
    // TODO: Do we want to go to full pivoting rather than partial pivoting for LU decomposition?
    // TODO: Consider checking Condition Number via LAPACK_zgecon
    // TODO: Consider ZGEEQU for equilibiration (scaling of rows / columns for better conditioning)

    // Calculate LU decomposition with partial pivoting
    triqs::arrays::vector<int> pivot_index_array(first_dim(g_mat_s1));
    int info = triqs::arrays::lapack::getrf(g_mat_s1, pivot_index_array, true);
    if (info != 0) {
      TRIQS_RUNTIME_ERROR << "lapack::getrf failed with code " << info;
    }

    // Extract det from LU decompositon (note permutations)
    dcomplex g_mat_s1_det = 1.0;
    int det_flip = 1;
    for (int i = 0; i < first_dim(g_mat_s1); i++) {
      g_mat_s1_det *= g_mat_s1(i, i);
      if (pivot_index_array(i) != i + 1) { // TODO: check +1 from fortran index
        det_flip = -det_flip;
      }
    }
    g_mat_s1_det *= det_flip;

    // Find cofactors by solving linear equation Ax = e1 * det(A)
    matrix<dcomplex> x_minors(first_dim(g_mat_s1), 1, FORTRAN_LAYOUT);
    x_minors() = 0;
    x_minors(0, 0) = g_mat_s1_det;

    info = triqs::arrays::lapack::getrs(g_mat_s1, x_minors, pivot_index_array);
    if (info != 0) {
      TRIQS_RUNTIME_ERROR << "lapack::getrs failed with code " << info;
    }

    // Multiply by Keldysh index parity & other spin determinant:
    x_minors *= GetBitParity(idx_kel) * determinant(g_mat_s2);

    // Bin: Find gf_index to bin to. We leave out first element, which connects external verticies only
    for (int i = 1; i < x_minors.size(); i++) {
      auto ind = all_config_1[col_pick_s1[i]];
      //result.accumulate(x_minors(i, 0), ind.contour.time, ind.contour.k_idx);
      result += x_minors(i, 0) * std::exp(1_j * omega * ind.contour.time);
    }
  }
  return std::make_pair(result, 1);
}

} // namespace keldy::impurity_oneband
