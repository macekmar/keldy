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

#include "wick_direct_time.hpp"
#include "keldy/common.hpp"
#include <complex>
#include <tuple>
#include <utility>

#pragma omp declare reduction(+ : keldy::dcomplex : omp_out += omp_in)

namespace {

inline int GetBit(int in, int offset) { return static_cast<int>((in & (1 << offset)) != 0); }

inline int GetBitParity(unsigned int in) { return 1 - 2 * __builtin_parity(in); }

} // namespace

namespace keldy::impurity_oneband {

// should we sort times?
std::pair<binner::sparse_binner_t<1>, int> integrand_g_direct_time::operator()(std::vector<double> const &times,
                                                                               bool const keep_u_hypercube) const {
  using namespace triqs::arrays;

  int order_n = times.size();
  if (order_n == 0) {
    TRIQS_RUNTIME_ERROR << "Order 0 is not implemented";
  }

  double time_binner = *std::min_element(times.cbegin(), times.cend());
  binner::sparse_binner_t<1> result;

  // Model is diagonal in spin
  if (external_A.spin != external_B.spin) {
    result.accumulate(0, time_binner);
    return std::make_pair(result, 0);
  }

  if (keep_u_hypercube) {
    // Integration starts a t = 0
    if (std::any_of(times.cbegin(), times.cend(), [](double t) { return t < 0.0; })) {
      result.accumulate(0, time_binner);
      return std::make_pair(result, 0);
    }
  }

  // copy for now since we change time-splitting
  auto a = external_A;
  auto b = external_B;
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

#pragma omp parallel for
  for (int i = 0; i < order_n; i++) {
    all_config_1[i] = gf_index_t{times[i], a.spin, forward, i};
    all_config_1[i + order_n] = gf_index_t{times[i], a.spin, backward, i};
    all_config_2[i] = gf_index_t{times[i], spin_t(1 - a.spin), forward, i};
    all_config_2[i + order_n] = gf_index_t{times[i], spin_t(1 - a.spin), backward, i};
  }

  // Index for external index in s1
  int external_idx = 2 * order_n;

  wick_matrix_s1(external_idx, external_idx) = g0(a, b, false);
#pragma omp parallel for
  for (int i = 0; i < 2 * order_n; i++) {
    wick_matrix_s1(external_idx, i) = g0(a, all_config_1[i]);
    wick_matrix_s1(i, external_idx) = g0(all_config_1[i], b);
    for (int j = 0; j < 2 * order_n; j++) {
      wick_matrix_s1(i, j) = g0(all_config_1[i], all_config_1[j]);
      wick_matrix_s2(i, j) = g0(all_config_2[i], all_config_2[j]);
    }
  }
  // std::cout << wick_matrix_s1 << std::endl;
  // std::cout << wick_matrix_s2 << std::endl;

  dcomplex integrand_result = 0.0;
  uint64_t nr_keldysh_configs = (uint64_t(1) << order_n);

  // Iterate over other Keldysh index configurations. Splict smaller determinant from precomuted matrix

#pragma omp parallel for reduction(+ : integrand_result)
  for (uint64_t idx_kel = 0; idx_kel < nr_keldysh_configs; idx_kel++) {
    // Indices of Rows / Cols to pick. Cycle through and shift by (0/1) * order_n depending on idx_kel configuration
    std::vector<int> col_pick_s2(order_n);
    for (int i = 0; i < order_n; i++) {
      col_pick_s2[i] = i + GetBit(idx_kel, i) * order_n;
    }
    std::vector<int> col_pick_s1 = col_pick_s2;
    col_pick_s1.push_back(external_idx);

    // Extract data into temporary matrices
    matrix<dcomplex> tmp_mat_s1(order_n + 1, order_n + 1);
    matrix<dcomplex> tmp_mat_s2(order_n, order_n);

    // std::cout << "col_pick_s1" << std::endl;

    // for(auto x: col_pick_s1) {
    //   std::cout << x << std::endl;
    // }

    // std::cout << "col_pick_s2" << std::endl;

    // for(auto x: col_pick_s2) {
    //   std::cout << x << std::endl;
    // }

    for (int i = 0; i < order_n + 1; ++i) {
      for (int j = 0; j < order_n + 1; ++j) {
        tmp_mat_s1(i, j) = wick_matrix_s1(col_pick_s1[i], col_pick_s1[j]);
      }
    }

    tmp_mat_s1(order_n, order_n) = 0; // cancelled by Keldysh indices

    for (int i = 0; i < order_n; ++i) {
      for (int j = 0; j < order_n; ++j) {
        tmp_mat_s2(i, j) = wick_matrix_s2(col_pick_s2[i], col_pick_s2[j]);
      }
    }
    integrand_result += GetBitParity(idx_kel) * determinant(tmp_mat_s1) * determinant(tmp_mat_s2);
  }

  // apply cutoff
  if (std::abs(integrand_result) < cutoff) {
    integrand_result = 0.;
  }

  // Multiply by overall factors (-1j) * (j)^n * (-1j)^n ?? FIXME
  result.accumulate(integrand_result, time_binner);
  return std::make_pair(result, 1);
}

} // namespace keldy::impurity_oneband
