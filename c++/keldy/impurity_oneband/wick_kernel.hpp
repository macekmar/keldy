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

#pragma once

#include "../common.hpp"
#include "model.hpp"

#include <mpi/mpi.hpp>
#include <mpi/vector.hpp>
#include <triqs/arrays/math_functions.hpp>
#include <triqs/arrays/mpi.hpp>
#include <triqs/utility/exceptions.hpp>
#include <triqs/utility/first_include.hpp>
#include <algorithm>
#include <vector>
#include <tuple>

namespace keldy::impurity_oneband {

using namespace triqs::arrays;
using namespace triqs::gfs;

// Small wrapper for data computed in kernel
class sparse_kernel_binner {
 public:
  std::vector<std::pair<gf_index_t, dcomplex>> data{};

  void bin_data(std::pair<gf_index_t, dcomplex> const & in){
    auto loc = std::find_if(std::begin(data), std::end(data), [&in](auto& el){return el.first == in.first;});
    if (loc != std::end(data)) {
        (*loc).second += in.second;
    } else {
      data.push_back(in);
    }
  }

  double sum_weights() {
    double result = 0.0;
    for (const auto &rh : data) {
      result += std::abs(rh.second);
    }
    return result;
  }

  template <typename T>
  sparse_kernel_binner &operator*=(T scalar) {
    for (auto &rh : data) {
      rh.second *= scalar;
    }
    return *this;
  }
};

/// Kernel binner for Green Function $K(Y, X')$ with binning happening over $Y$
/// and $X'$ fixed by boundary conditions.
///
/// TODO: How to include spin up / down separatley.
class kernel_binner {
  double t_min = 0.0;
  double t_max = 1.0;
  int n_bins = 100;
  double bin_size = 0.01;

  int nr_point_dropped = 0;

  array<dcomplex, 2> values;    // bin (time), keldysh index, orbitals
  array<uint64_t, 2> nr_values; // bin (time), keldysh index, orbitals
  array<double, 1> bin_times;

 public:
  kernel_binner() : kernel_binner(0.0, 1.0, 100) {};
  kernel_binner(double t_min_, double t_max_, int n_bins_)
     : t_min(t_min_),
       t_max(t_max_),
       n_bins(n_bins_),
       //  external_point_X(std::move(external_point_X_)),
       values(n_bins_, 2),
       nr_values(n_bins_, 2),
       bin_times(n_bins_) {
    values() = 0;
    nr_values() = 0;

    // make checks?
    bin_size = (t_max - t_min) / n_bins;

    for (int i = 0; i < n_bins; ++i) {
      bin_times(i) = t_min + (i + 0.5) * bin_size;
    }
  }

  auto get_values() const { return values; }
  auto get_nr_values() const { return nr_values; }
  auto get_bin_times() const { return bin_times; }
  auto get_bin_size() const { return bin_size; }
  auto get_nr_point_dropped() const { return nr_point_dropped; }

  /// Includes boundary points, so t_min <= t <= t_max. t_max gets put in last bin
  // CHECK THIS TO BE CORRECT!
  void accumulate(gf_index_t const &a, dcomplex value) {
    if (t_min <= a.time && a.time < t_max) {
      int bin = int((a.time - t_min) / bin_size);
      values(bin, a.k_idx) += value;
      nr_values(bin, a.k_idx)++;
    } else if (a.time == t_max) {
      values(n_bins - 1, a.k_idx) += value;
      nr_values(n_bins - 1, a.k_idx)++;
    } else {
      nr_point_dropped++;
    }
  }

  kernel_binner &operator+=(sparse_kernel_binner const &rhs) {
    for (const auto &rh : rhs.data) {
      accumulate(rh.first, rh.second);
    }
    return *this;
  }

  template <typename T>
  kernel_binner &operator*=(T scalar) {
    values *= scalar;
    return *this;
  }

  friend kernel_binner mpi_reduce(kernel_binner const &a, mpi::communicator c, int root, bool all, MPI_Op op);
};



inline kernel_binner mpi_reduce(kernel_binner const &in, mpi::communicator c = {}, int root = 0, bool all = false,
                         MPI_Op op = MPI_SUM) {
  if(op != MPI_SUM){
    TRIQS_RUNTIME_ERROR << "mpi_reduce of kernel_binner can only be performed with op = MPI_SUM";
  }

  auto t_min_vec = mpi::gather(std::vector<double>({in.t_min}), c, root, all);
  auto t_max_vec = mpi::gather(std::vector<double>({in.t_max}), c, root, all);
  auto n_bins_vec = mpi::gather(std::vector<int>({in.n_bins}), c, root, all);

  bool all_eq = std::equal(t_min_vec.begin() + 1, t_min_vec.end(), t_min_vec.begin());
  all_eq = all_eq && std::equal(t_max_vec.begin() + 1, t_max_vec.end(), t_max_vec.begin());
  all_eq = all_eq && std::equal(n_bins_vec.begin() + 1, n_bins_vec.end(), n_bins_vec.begin());

  if(!all_eq){
    TRIQS_RUNTIME_ERROR << "kernel_binner to mpi_reduce are not defined with same bins. This is not supported!";
  }

  kernel_binner out;
  // since bins are the same, just copy local values
  out.t_min = in.t_min;
  out.t_max = in.t_max;
  out.n_bins = in.n_bins;
  out.bin_size = in.bin_size;
  out.bin_times = in.bin_times;
  out.nr_point_dropped = 0;

  if (!all) {
    out.values = mpi::reduce(in.values, c, 0, false, MPI_SUM);
    out.nr_values = mpi::reduce(in.nr_values, c, 0, false, MPI_SUM);
    out.nr_point_dropped = mpi::reduce(in.nr_point_dropped, c, 0, false, MPI_SUM);
  } else {
    out.values = mpi::reduce(in.values, c, 0, true, MPI_SUM);
    out.nr_values = mpi::reduce(in.nr_values, c, 0, true, MPI_SUM);
    out.nr_point_dropped = mpi::reduce(in.nr_point_dropped, c, 0, true, MPI_SUM);
  }
  return out;
}


// template<kernel_binner>
// constexpr bool is_binned_variable = true;

class integrand_g_kernel {
  g0_keldysh_contour_t g0;
  gf_index_t g_idx_X; // Fixed Point in Kernal

  // bool expand_col = true; // expand_row = false
  // double condition_numebr_tol;

 public:
  /// Returns integrand for the specified times
  using result_t = sparse_kernel_binner;
  result_t operator()(std::vector<double> const &times) const;

  integrand_g_kernel(g0_keldysh_contour_t g0_, gf_index_t g_idx_X_)
     : g0(std::move(g0_)), g_idx_X(std::move(g_idx_X_)){};
};

} // namespace keldy::impurity_oneband
