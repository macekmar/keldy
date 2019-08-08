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

namespace {

inline int GetBit(int in, int offset) { return (in & (1 << offset)) != 0; }

inline int GetBitParity(unsigned int in) { return 1 - 2 * __builtin_parity(in); }

} // namespace

namespace keldy::impurity_oneband {

// should we sort times?
dcomplex integrand_g_t1t2_direct::operator()(std::vector<double> const &times) const {
  using namespace triqs::arrays;
  int order_n = times.size();

  auto &a = external_A;
  auto &b = external_B;

  if (a.spin != b.spin) {
    return 0.0;
  }

  if (order_n == 0) {
    return g0(a, b, order_n, order_n);
  }

  if(*std::min_element(times.begin(), times.end()) < 0){ // can replace with a for_any
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
    current_config_1[i] = gf_index_t{times[i], a.spin, forward};
    current_config_2[i] = gf_index_t{times[i], spin_t(1 - a.spin), forward};
  }

  wick_matrix_s1(order_n, order_n) = g0(a, b, order_n, order_n); // equal if we set == 0;
  for (int i = 0; i < order_n; i++) {
    wick_matrix_s1(order_n, i) = g0(a, current_config_1[i], order_n, i);
    wick_matrix_s1(i, order_n) = g0(current_config_1[i], b, i, order_n);

    for (int j = 0; j < order_n; j++) {
      wick_matrix_s1(i, j) = g0(current_config_1[i], current_config_1[j], i, j);
      wick_matrix_s2(i, j) = g0(current_config_2[i], current_config_2[j], i, j);
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
    wick_matrix_s1(order_n, nlc) = g0(a, current_config_1[nlc], order_n, nlc);
    wick_matrix_s1(nlc, order_n) = g0(current_config_1[nlc], b, nlc, order_n);

    for (int i = 0; i < order_n; i++) {
      wick_matrix_s1(i, nlc) = g0(current_config_1[i], current_config_1[nlc], i, nlc);
      wick_matrix_s1(nlc, i) = g0(current_config_1[nlc], current_config_1[i], nlc, i);

      wick_matrix_s2(i, nlc) = g0(current_config_2[i], current_config_2[nlc], i, nlc);
      wick_matrix_s2(nlc, i) = g0(current_config_2[nlc], current_config_2[i], nlc, i);
    }

    integrand_result +=
       parity * triqs::arrays::determinant(wick_matrix_s1) * triqs::arrays::determinant(wick_matrix_s2);
  }
  return integrand_result;
}

} // namespace keldy::impurity_oneband
