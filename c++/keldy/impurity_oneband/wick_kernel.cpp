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

#include "wick_kernel.hpp"

namespace {

inline int GetBit(int in, int offset) { return (in & (1 << offset)) != 0; }

inline int GetBitParity(unsigned int in) { return 1 - 2 * __builtin_parity(in); }

} // namespace

namespace keldy::impurity_oneband {

// TODO: Where is sorting of times done?
// TODO: Efficnencies for spin symmetric case.

void integrand_g_kernel::operator()(kernel_binner &kernel, std::vector<double> const &times) const {
  using namespace triqs::arrays;
  // Interaction starts a t = 0
  if (*std::min_element(times.begin(), times.end()) < 0) { // can replace with a for_any
    return;
  }

  int order_n = times.size();

  // if (order_n == 0) {
  //   return 1.0; // TODO: CHECK NORMALIZATION
  // }

  auto a = external_point_X; // This is now a dummy index
  auto b = external_point_X;
  // define time-splitting for external-points
  a.timesplit_n = order_n;
  b.timesplit_n = order_n;

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

  dcomplex integrand_result = 0.0;
  uint64_t nr_keldysh_configs = (uint64_t(1) << order_n);

  // Iterate over other Keldysh index configurations. Splict smaller determinant from precomuted matrix
  // #pragma omp parallel for reduce(+: integrand_result)
  for (uint64_t idx_kel = 0; idx_kel < nr_keldysh_configs - 1; idx_kel++) {
    // Indices of Rows / Cols to pick. Cycle through and shift by (0/1) * order_n depending on idx_kel configuration
    std::vector<int> col_pick_s1(order_n + 1);
    std::vector<int> col_pick_s2(order_n);

    col_pick_s1[0] = external_idx;
    for (int i = 0; i < order_n; i++) {
      col_pick_s1[i + 1] = = i + GetBit(idx_kel, i) * order_n;
      col_pick_s2[i] = col_pick_s1[i + 1];
    }

    // Extract data into temporar matrices
    matrix<dcomplex> tmp_mat_s1(order_n + 1, order_n + 1);
    matrix<dcomplex> tmp_mat_s2(order_n, order_n);

    for (auto [i, j] : itertools::zip(col_pick_s2, col_pick_s2)) {
      tmp_mat_s2(i, j) = wick_matrix_s2(col_pick_s2[i], col_pick_s2[j]);
    }

    for (auto [i, j] : itertools::zip(col_pick_s1, col_pick_s1)) {
      tmp_mat_s1(i, j) = wick_matrix_s1(col_pick_s1[i], col_pick_s1[j]);
    }

    array<int, 1> pivot_index_array(tmp_mat_s1.rows());
    pivot_index_array() = 0;

    auto tmp_mat_s1_copy = tmp_mat_s1;
    int info = 0;

    // TODO: FIX LAPACK WRAPPING / SIGNATURE:
    // TODO: Do we want to go to full pivoting rather than partial pivoting?

    // Calculate LU decomposition with partial pivoting
    LAPACK_zgetrf(&tmp_mat_s1.rows(), &tmp_mat_s1.cols(), tmp_mat_s1.data(), &tmp_mat_s1.rows(), tmp_mat_s1.data(),
                  &info);
    if (info != 0) TRIQS_RUNTIME_ERROR << "Lapack zgetrf failed with error code " << info;

    // TODO: Calculate Condition Number:
    // LAPACK_zgecon(char const* norm, lapack_int const* n,lapack_complex_double const* A, lapack_int const* lda,
    // double const* anorm,ouble* rcond,lapack_complex_double* work,double* rwork,lapack_int* info )

    // TODO: ZGEEQU for equilibiration (scaling of rows / columns for better conditioning)

    // Extract det from LU decompositon (note permutations)
    dcomplex tmp_mat_s1_det = 1.0;
    for (int i = 0; i < tmp_mat_s1.rows(); i++) {
      tmp_mat_s1_det *= tmp_mat_s1(i, i);
      if (pivot_index_array(i) != i) {
        tmp_mat_s1_det *= 1;
      }
    }

    // Find cofactors by solving linear equation Ax = e1 * det(A)
    array<dcomplex, 1> x_minors(tmp_mat_s1.rows());
    x_minors() = 0;
    x_minors(0) = tmp_mat_s1_det;

    int nrhs = 1;
    LAPACK_zgetrs('N', &tmp_mat_s1.rows(), &nrhs, tmp_mat_s1.data(), &tmp_mat_s1.rows(), tmp_mat_s1.data(),
                  x_minors.data(), &x_minors.rows(), &info);

    // Multiply by signs / parity / other spin determinant:
    x_minors *= GetBitParity(idx_kel) * determinant(tmp_mat_s2);

    // Bin: Find gf_index to bin to. We leave out first element, which connects external verticies only
    for (int i = 1; i < x_minors.size(); i++) {
      kernel.accumulate(all_config_1[col_pick_s1[i]], x_minors(i));
    }
  }
}

} // namespace keldy::impurity_oneband
