/*******************************************************************************
 *
 * keldy
 *
 * Copyright (c) 2020 The Simons Foundation
 * Copyright (c) 2020 CEA: Commissariat à l’énergie atomique
 *                         et aux énergies alternatives
 *   authors: Philipp Dumitrescu, Marjan Macek, Corentin Bertrand
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

#include "wick_direct.hpp"
#include <complex>
#include <tuple>
#include <utility>

#pragma omp declare reduction(+ : dcomplex : omp_out += omp_in)

namespace {

inline int GetBit(int in, int offset) { return static_cast<int>((in & (1 << offset)) != 0); }

inline int GetBitParity(unsigned int in) { return 1 - 2 * __builtin_parity(in); }

} // namespace

namespace keldy::impurity_oneband {



void green_function_config (std::vector<double> const &times, spin_t const spin, std::vector<gf_index_t> &config) {

  int order = times.size();
  assert (config.size() == 2 * order);

  for (int i = 0; i < order; i++) {
    config[i] = gf_index_t{times[i], spin, forward, i};
    config[i + order] = gf_index_t{times[i], spin, backward, i};
  }
}


template <typename T>
void reorder_matrix (T const &mat_in, T &mat_out, std::vector<int> const &ordering) {

  int size = ordering.size();
  // check matrix size;

  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      mat_out(i, j) = mat_in(ordering[i], ordering[j]);
    }
  }
}


template <typename T>
void calc_wick_matrix (g0_keldysh_contour_t g0, std::vector<gf_index_t> const &config, T &wickmat, int ext_idx=-1, gf_index_t a=gf_index_t(), gf_index_t b=gf_index_t()) {

  int size = config.size();

  if (ext_idx != -1) {
    wickmat(ext_idx, ext_idx) = g0(a, b, false);
    #pragma omp parallel for
    for (int i = 0; i < size; i++) {
      wickmat(ext_idx, i) = g0(a, config[i]);
      wickmat(i, ext_idx) = g0(config[i], b);
    }
  }
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      wickmat(i, j) = g0(config[i], config[j]);
    }
  }
}



void calc_ordering_vector(int const keldysh_idx, int const order, std::vector<int> &ordering){
  assert (ordering.size() == order);
  for (int i = 0; i < order; i++) {
    ordering[i] = i + GetBit(keldysh_idx, i) * order;
  }
}



template <typename T>
class WickMatrix {

  std::vector<double> times;
  int order;
  gf_index_t a, b;
  g0_keldysh_contour_t g0;

 public:

  uint64_t max_keldysh_configs;

  WickMatrix(std::vector<double> const & _times, g0_keldysh_contour_t _g0, gf_index_t _a, gf_index_t _b) : 
  times{_times}, g0{_g0}, a{_a}, b{_b} {

  order = times.size();
  assert (order != 0);
  size = 2 * order;

  max_keldysh_configs = (uint64_t(1) << order);
  external_idx = size;

  config1.resize(size);
  config2.resize(size);
  ordering_1.resize(size + 1);
  ordering_2.resize(size);

    for (int i = 0; i < order; i++) {
      ordering_2[i] = 42;
    }


  wick_mat_1.resize(size + 1, size + 1);
  wick_mat_2.resize(size, size);
  wick_mat_reord_1.resize(size + 1, size + 1);
  wick_mat_reord_2.resize(size, size);

  green_function_config(times, a.spin, config1);
  green_function_config(times, spin_t(1 - a.spin), config2);

  calc_wick_matrix<T>(g0, config1, wick_mat_1, external_idx, a, b);
  calc_wick_matrix<T>(g0, config2, wick_mat_2);

  }

  dcomplex determinant(int keldysh_idx){
  // Calculate the determinant of the Wick matrix.
    assert (keldysh_idx >= _0);
    assert (keldysh_idx < max_keldysh_configs);

    permutation(keldysh_idx);
    return GetBitParity(keldysh_idx) * triqs::arrays::determinant(wick_mat_reord_1) * triqs::arrays::determinant(wick_mat_reord_2);
  }

  dcomplex permanent(int keldysh_idx){
  // Calculate the permanent of the Wick matrix.
  return 0.;
  }

  dcomplex kernel(int keldysh_idx){
  // Calculate the kernel of the Wick matrix.
  return 0.;
  }

  std::vector<int> get_ordering2(int keldysh_idx){

    assert (keldysh_idx >= _0);
    assert (keldysh_idx < max_keldysh_configs);

    std::vector<int> xx(order);

    for (int i = 0; i < order; i++) {
        xx[i] = i + GetBit(keldysh_idx, i) * order;
     };
    return xx;
  }

  void permutation(int keldysh_idx){
  // Perform a permutation of the Wick matrix.
    assert (keldysh_idx >= _0);
    assert (keldysh_idx < max_keldysh_configs);

    std::vector<int> xx(order);


    //calc_ordering_vector(keldysh_idx, order, ordering_2);

    for (int i = 0; i < order; i++) {
       // std::cout << "inside i: " << i << std::endl;
       // std::cout << "inside routine: " << ordering_2[i] << "  " << i + GetBit(keldysh_idx, i) * order << std::endl;
        xx[i] = i + GetBit(keldysh_idx, i) * order;
        ordering_2[i] = xx[i];
     };

      std::cout << "inside routine: " << ordering_2[0] << std::endl;
      std::cout << std::flush;

//    for (int i = 0; i < order; i++) {
//      ordering_1[i] = ordering_2[i];
//    }
//    ordering_1[order + 1] = external_idx;

//    reorder_matrix<T> (wick_mat_1, wick_mat_reord_1, ordering_1);
//    reorder_matrix<T> (wick_mat_2, wick_mat_reord_2, ordering_2);
  }


  int external_idx;
  int size;
  std::vector<gf_index_t> config1;
  std::vector<gf_index_t> config2;
  std::vector<int> ordering_1;
  std::vector<int> ordering_2;
  T wick_mat_1;
  T wick_mat_2;
  T wick_mat_reord_1;
  T wick_mat_reord_2;

};

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> &input)
{
	for (auto const& i: input) {
		os << i << " ";
	}
	return os;
}


// should we sort times?
std::pair<dcomplex, int> integrand_g_direct::operator()(std::vector<double> const &times,
                                                        bool const keep_u_hypercube) const {
  using namespace triqs::arrays;

  // Model is diagonal in spin
  if (external_A.spin != external_B.spin) {
    return std::make_pair(0.0, 0);
  }

  // Interaction starts a t = 0
  if (keep_u_hypercube) {
    if (std::any_of(times.cbegin(), times.cend(), [](double t) { return t < 0.0; })) {
      return std::make_pair(0.0, 0);
    }
  }

  int order_n = times.size();

  // copy for now since we change time-splitting
  auto a = external_A;
  auto b = external_B;
  // define time-splitting for external-points
  a.contour.timesplit_n = order_n;
  b.contour.timesplit_n = order_n;

  if (order_n == 0) {
    return std::make_pair(g0(a, b, false), 1);
  }

  auto wick_matrix = WickMatrix<triqs::arrays::matrix<dcomplex>> (times, g0, a, b);


  // Pre-Comute Large Matrix.
  // "s1": Same spin as external indices / "s2": Opposite spin
//  matrix<dcomplex> wick_matrix_s1(2 * order_n + 1, 2 * order_n + 1);
//  matrix<dcomplex> wick_matrix_s2(2 * order_n, 2 * order_n);

  // Vector of indices for Green functions
//  std::vector<gf_index_t> all_config_1(2 * order_n);
//  std::vector<gf_index_t> all_config_2(2 * order_n);

//  green_function_config(times, a.spin, all_config_1);
//  green_function_config(times, spin_t(1 - a.spin), all_config_2);

  // Index for external index in s1
  int external_idx = 2 * order_n;

//  calc_wick_matrix<matrix<dcomplex>>(g0, all_config_1, wick_matrix_s1, external_idx, a, b);
//  calc_wick_matrix<matrix<dcomplex>>(g0, all_config_2, wick_matrix_s2);

  auto wick_matrix_s1 = wick_matrix.wick_mat_1;
  auto wick_matrix_s2 = wick_matrix.wick_mat_2;

  dcomplex integrand_result = 0.0;
  uint64_t nr_keldysh_configs = (uint64_t(1) << order_n);

  // Iterate over other Keldysh index configurations. Splict smaller determinant from precomuted matrix

#pragma omp parallel for reduction(+ : integrand_result)
  for (uint64_t idx_kel = 0; idx_kel < wick_matrix.max_keldysh_configs; idx_kel++) {
    // Indices of Rows / Cols to pick. Cycle through and shift by (0/1) * order_n depending on idx_kel configuration

    std::vector<int> col_pick_s2(order_n);
    calc_ordering_vector(idx_kel, order_n, col_pick_s2);
    std::vector<int> col_pick_s1 = col_pick_s2;
    col_pick_s1.push_back(external_idx);

//    col_pick_s1 = wick_matrix.ordering_1;
    wick_matrix.permutation(idx_kel);
    //auto order2 = wick_matrix.ordering_2;



    int i = 1;

        wick_matrix.permutation(idx_kel);
          auto elem = wick_matrix.ordering_2[i];
        wick_matrix.permutation(idx_kel);
          auto elem1 = wick_matrix.ordering_2;
        wick_matrix.permutation(idx_kel);
          auto elem2 = wick_matrix.ordering_2;
        wick_matrix.permutation(idx_kel);
          auto elem3 = wick_matrix.ordering_2;
        wick_matrix.permutation(idx_kel);
          auto elem4 = wick_matrix.ordering_2;

          auto elem5 = wick_matrix.get_ordering2(idx_kel);

      if (elem != col_pick_s2[i]){
        std::cout << "error detected " << wick_matrix.ordering_2[i] << " " << col_pick_s2[i] <<std::endl;
      std::cout << std::flush;
        //print<int>(col_pick_s2);
        std::cout << "elem1: " << elem1 << std::endl;
        std::cout << "elem2: " << elem2 << std::endl;
        std::cout << "elem3: " << elem3 << std::endl;
        std::cout << "elem4: " << elem4 << std::endl;
        //std::cout << "elem5: " << elem5 << std::endl;
        std::cout << "i: " << i << std::endl;
        std::cout << "kel_idx: " << idx_kel << std::endl;
      std::cout << std::flush;
//        wick_matrix.permutation(idx_kel);
        std::cout << "is: " << wick_matrix.ordering_2[i] << std::endl;
        std::cout << "ref: " << col_pick_s2[i] << std::endl;
        std::cout << "ref2: " <<i + GetBit(idx_kel, i) * order_n << std::endl;
        std::cout << "error abort " << wick_matrix.ordering_2[i] << " " << col_pick_s2[i] << std::endl;
        std::cout << std::flush;
        abort();
        //};
      //std::cout << "i: " << i << " ismref: " << order2[i] - col_pick_s2[i] << std::endl;
    }

    // Extract data into temporary matrices
    matrix<dcomplex> tmp_mat_s1(order_n + 1, order_n + 1);
    matrix<dcomplex> tmp_mat_s2(order_n, order_n);


    reorder_matrix<matrix<dcomplex>> (wick_matrix_s1, tmp_mat_s1, col_pick_s1);
    reorder_matrix<matrix<dcomplex>> (wick_matrix_s2, tmp_mat_s2, col_pick_s2);

    integrand_result += GetBitParity(idx_kel) * determinant(tmp_mat_s1) * determinant(tmp_mat_s2);
    //integrand_result += wick_matrix.determinant(idx_kel);

  }

  // apply cutoff
  if (std::abs(integrand_result) < cutoff) {
    integrand_result = 0.;
  }

  // Multiply by overall factors (-1j) * (j)^n * (-1j)^n ?? FIXME
  return std::make_pair(integrand_result, 1);
}

// *******************************************************

// Old Method Based on Sequential Constructions
// Copy a,b as need to modify time-splitting
dcomplex integrand_g_direct_grey(gf_index_t a, gf_index_t b, g0_keldysh_contour_t const &g0,
                                 std::vector<double> const &times) {

  // TODO: should we sort times?

  using namespace triqs::arrays;
  int order_n = times.size();

  // copy for now since we change time-splitting

  if (a.spin != b.spin) {
    return 0.0;
  }

  // define time-splitting for external-points
  a.contour.timesplit_n = order_n;
  b.contour.timesplit_n = order_n;

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
    current_config_1[nlc].contour.k_idx = keldysh_idx_t(1 - current_config_1[nlc].contour.k_idx);
    current_config_2[nlc].contour.k_idx = keldysh_idx_t(1 - current_config_2[nlc].contour.k_idx);

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
