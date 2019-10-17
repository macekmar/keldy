// //******************************************************************************
// //
// // keldy
// //
// // Copyright (C) 2019, The Simons Foundation
// // authors: Philipp Dumitrescu
// //
// // keldy is free software: you can redistribute it and/or modify it under the
// // terms of the GNU General Public License as published by the Free Software
// // Foundation, either version 3 of the License, or (at your option) any later
// // version.
// //
// // keldy is distributed in the hope that it will be useful, but WITHOUT ANY
// // WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// // FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// // details.
// //
// // You should have received a copy of the GNU General Public License along with
// // keldy. If not, see <http://www.gnu.org/licenses/>.
// //
// //******************************************************************************

// #pragma once

// #include "../common.hpp"
// #include "model.hpp"

// #include <triqs/utility/first_include.hpp>
// // #include <triqs/utility/variant.hpp>
// #include <vector>

// namespace keldy::impurity_oneband {

// using namespace triqs::arrays;
// using namespace triqs::gfs;

// // Small wrapper for data
// class sparse_kernel {
// };

// /// Kernel binner for Green Function $K(Y, X')$ with binning happening over $Y$
// /// and $X'$ fixed by boundary conditions.
// ///
// /// TODO: How to include spin up / down separatley.
// ///
// class kernel_binner {
//   double t_min = 0.0;
//   double t_max = 1.0;
//   int nr_bins = 100;
//   double bin_size = 0.01;

//   // // External point $X'$
//   // gf_index_t external_point_X;

//   int nr_point_dropped = 0;

//   array<dcomplex, 2> values;    // bin (time), keldysh index, orbitals
//   array<uint64_t, 2> nr_values; // bin (time), keldysh index, orbitals
//   array<double, 1> bin_times;

//  public:
//   kernel_binner() = default;
//   kernel_binner(double t_min_, double t_max_, int nr_bins_, gf_index_t external_point_X_)
//      : t_min(t_min_),
//        t_max(t_max_),
//        nr_bins(nr_bins_),
//       //  external_point_X(std::move(external_point_X_)),
//        values(nr_bins, 2),
//        nr_values(nr_bins, 2),
//        bin_times(nr_bins) {

//     values() = 0;
//     nr_values() = 0;

//     // make checks?
//     bin_size = (t_max - t_min) / nr_bins;

//     for (int i = 0; i < nr_bins; ++i) {
//       bin_times(i) = t_min + (i + 0.5) * bin_size;
//     }
//   }

//   auto get_values() const { return values; }   
//   auto get_nr_values() const { return nr_values; }
//   auto get_bin_times() const { return bin_times; }
//   auto get_bin_size() const { return bin_size; }
//   auto get_nr_point_dropped() const { return nr_point_dropped; }

//   /// Includes boundary points, so t_min <= t <= t_max. t_max gets put in last bin
//   void accumulate(gf_index_t const &a, dcomplex value) {
//     if (t_min < a.time && a.time < t_max) {
//       int bin = int((a.time - t_min) / bin_size);
//       values(bin, a.k_idx) += value;
//       nr_values(bin, a.k_idx)++;
//     } else if (a.time == t_max) {
//       values(nr_bins - 1, a.k_idx) += value;
//       nr_values(nr_bins - 1, a.k_idx)++;
//     } else {
//       nr_point_dropped++;
//     }
//   }

//   kernel_binner& operator+=(std::vector<std::pair<gf_index_t, dcomplex>> const & rhs){
//     for(int i = 0; i < rhs.size(); i++){
//       accumulate(rhs[i].first, rhs[i].second);
//     }
//     return *this;
//   }

// };

// // implement this: double * std::vector<std::pair<gf_index_t, dcomplex>>



// class integrand_g_kernel {
//   g0_keldysh_contour_t g0;
//   gf_index_t g_idx_X; // Fixed Point in Kernal

//   // bool expand_col = true; // expand_row = false
//   // double condition_numebr_tol;

//  public:
//   /// Returns integrand for the specified times
//   std::vector<std::pair<gf_index_t, dcomplex>> operator()(std::vector<double> const &times) const;

//   integrand_g_kernel(g0_keldysh_contour_t g0_, gf_index_t g_idx_X_) : g0(std::move(g0_)), g_idx_X(std::move(g_idx_X_)) {};
// };

// // template<kernel_binner>
// // constexpr bool is_binned_variable = true;


// } // namespace keldy::impurity_oneband

