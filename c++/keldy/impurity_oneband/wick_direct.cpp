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

#include "wick_direct.hpp"
#include <complex>

#pragma omp declare reduction(+ : dcomplex : omp_out += omp_in)

namespace {

inline int GetBit(int in, int offset) { return static_cast<int>((in & (1 << offset)) != 0); }

inline int GetBitParity(unsigned int in) { return 1 - 2 * __builtin_parity(in); }

} // namespace

namespace keldy::impurity_oneband {

// should we sort times?
dcomplex integrand_g_t1t2_direct::operator()(std::vector<double> const &times) const {
  using namespace triqs::arrays;

  // Model is diagonal in spin
  if (external_A.spin != external_B.spin) {
    return 0.0;
  }

  // Interaction starts a t = 0
  if (*std::min_element(times.begin(), times.end()) < 0) { // can replace with a for_any
    return 0.0;
  }

  int order_n = times.size();

  // copy for now since we change time-splitting
  auto a = external_A;
  auto b = external_B;
  // define time-splitting for external-points
  a.timesplit_n = order_n;
  b.timesplit_n = order_n;

  if (order_n == 0) {
    return g0(a, b, false);
  }

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

  #pragma omp parallel for reduction(+:integrand_result)
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
    for (int i = 0; i < order_n; ++i) {
      for (int j = 0; j < order_n; ++j) {
        tmp_mat_s2(i, j) = wick_matrix_s2(col_pick_s2[i], col_pick_s2[j]);
      }
    }
    integrand_result += GetBitParity(idx_kel) * determinant(tmp_mat_s1) * determinant(tmp_mat_s2);
  }

  // Multiply by overall factors (-1j) * (j)^n * (-1j)^n ?? FIXME
  return integrand_result;
}

// *******************************************************

// Old Method Based on Sequential Constructions
// Copy a,b as need to modify time-splitting
dcomplex integrand_g_t1t2_direct_grey(gf_index_t a, gf_index_t b, g0_keldysh_contour_t const &g0,
                                      std::vector<double> const &times) {

  // TODO: should we sort times?

  using namespace triqs::arrays;
  int order_n = times.size();

  // copy for now since we change time-splitting

  if (a.spin != b.spin) {
    return 0.0;
  }

  // define time-splitting for external-points
  a.timesplit_n = order_n;
  b.timesplit_n = order_n;

  if (order_n == 0) {
    return g0(a, b, false);
  }

  if (*std::min_element(times.begin(), times.end()) < 0) { // can replace with a for_any
    return 0.0;
  }

  // must allow to be flippable
  matrix<dcomplex> wick_matrix_s1(order_n + 1, order_n + 1);
  matrix<dcomplex> wick_matrix_s2(order_n, order_n);

  // Construct for initial configuration (all forward contour):

  // Vector of gf vertex
  std::vector<gf_index_t> current_config_1(order_n);
  std::vector<gf_index_t> current_config_2(order_n);
  for (int i = 0; i < order_n; i++) {
    current_config_1[i] = gf_index_t{times[i], a.spin, forward, i};
    current_config_2[i] = gf_index_t{times[i], spin_t(1 - a.spin), forward, i};
  }

  wick_matrix_s1(order_n, order_n) = g0(a, b, false);
  for (int i = 0; i < order_n; i++) {
    wick_matrix_s1(order_n, i) = g0(a, current_config_1[i]);
    wick_matrix_s1(i, order_n) = g0(current_config_1[i], b);

    for (int j = 0; j < order_n; j++) {
      wick_matrix_s1(i, j) = g0(current_config_1[i], current_config_1[j]);
      wick_matrix_s2(i, j) = g0(current_config_2[i], current_config_2[j]);
    }
  }

  dcomplex integrand_result = determinant(wick_matrix_s1) * determinant(wick_matrix_s2);
  // TRIQS_PRINT(integrand_result);

  int parity = 1;
  uint64_t nr_keldysh_configs = (uint64_t(1) << order_n);
  // Iterate over other Keldysh index configurations
  for (uint64_t idx_kel = 0; idx_kel < nr_keldysh_configs - 1; idx_kel++) {
    // Use Grey code to select a single row / column to flip [cf Numerical Recipes, Sect.20.2]
    // ffs starts at 1, returns the position of the 1st (least significant) bit set to 1. ~n has bites inversed compared with n.
    int nlc = (idx_kel < nr_keldysh_configs - 1 ? ffs(~idx_kel) : order_n) - 1;

    parity = -parity;

    // implicit cast
    current_config_1[nlc].k_idx = keldysh_idx_t(1 - current_config_1[nlc].k_idx);
    current_config_2[nlc].k_idx = keldysh_idx_t(1 - current_config_2[nlc].k_idx);

    // Connect to External Vertices
    wick_matrix_s1(order_n, nlc) = g0(a, current_config_1[nlc]);
    wick_matrix_s1(nlc, order_n) = g0(current_config_1[nlc], b);

    for (int i = 0; i < order_n; i++) {
      wick_matrix_s1(i, nlc) = g0(current_config_1[i], current_config_1[nlc]);
      wick_matrix_s1(nlc, i) = g0(current_config_1[nlc], current_config_1[i]);

      wick_matrix_s2(i, nlc) = g0(current_config_2[i], current_config_2[nlc]);
      wick_matrix_s2(nlc, i) = g0(current_config_2[nlc], current_config_2[i]);
    }

    integrand_result += parity * determinant(wick_matrix_s1) * determinant(wick_matrix_s2);
  }
  return integrand_result;
}

} // namespace keldy::impurity_oneband
