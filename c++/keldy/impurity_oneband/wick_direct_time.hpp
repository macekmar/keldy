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

#include <triqs/utility/first_include.hpp>
#include <triqs/utility/variant.hpp>
#include <vector>

namespace keldy::impurity_oneband {

struct singleton_binner {
 public:
  double time;
  dcomplex value;

  template <typename T>
  singleton_binner &operator*=(T rhs) {
    this->value *= rhs;
    return *this;
  };
};

class binner_1d {
  double t_min = 0.0;
  double t_max = 1.0;
  int n_bins = 100;
  double bin_size = 0.01;

  int nr_point_dropped = 0;

  array<dcomplex, 1> values;
  array<uint64_t, 1> nr_values;
  array<double, 1> bin_times;

 public:
  binner_1d() : binner_1d(0.0, 1.0, 100){};
  binner_1d(double t_min_, double t_max_, int n_bins_)
     : t_min(t_min_), t_max(t_max_), n_bins(n_bins_), values(n_bins_), nr_values(n_bins_), bin_times(n_bins_) {
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
  void accumulate(double time, dcomplex value) {
    if (t_min <= time && time < t_max) {
      int bin = int((time - t_min) / bin_size);
      values(bin) += value;
      nr_values(bin)++;
    } else if (time == t_max) {
      values(n_bins - 1) += value;
      nr_values(n_bins - 1)++;
    } else {
      nr_point_dropped++;
    }
  }

  binner_1d &operator+=(singleton_binner const &rhs) {
    accumulate(rhs.time, rhs.value);
    return *this;
  }

  template <typename T>
  binner_1d &operator*=(T scalar) {
    values *= scalar;
    return *this;
  }

  friend binner_1d mpi_reduce(binner_1d const &a, mpi::communicator c, int root, bool all, MPI_Op op);
};

inline binner_1d CPP2PY_IGNORE mpi_reduce(binner_1d const &in, mpi::communicator c = {}, int root = 0, bool all = false,
                                          MPI_Op op = MPI_SUM) {
  if (op != MPI_SUM) {
    TRIQS_RUNTIME_ERROR << "mpi_reduce of binner_1d can only be performed with op = MPI_SUM";
  }

  auto t_min_vec = mpi::gather(std::vector<double>({in.t_min}), c, root, all);
  auto t_max_vec = mpi::gather(std::vector<double>({in.t_max}), c, root, all);
  auto n_bins_vec = mpi::gather(std::vector<int>({in.n_bins}), c, root, all);

  bool all_eq = std::equal(t_min_vec.begin() + 1, t_min_vec.end(), t_min_vec.begin());
  all_eq = all_eq && std::equal(t_max_vec.begin() + 1, t_max_vec.end(), t_max_vec.begin());
  all_eq = all_eq && std::equal(n_bins_vec.begin() + 1, n_bins_vec.end(), n_bins_vec.begin());

  if (!all_eq) {
    TRIQS_RUNTIME_ERROR << "binner_1d to mpi_reduce are not defined with same bins. This is not supported!";
  }

  binner_1d out;
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

//-----------------------

class integrand_g_direct_time {
  g0_keldysh_contour_t g0;
  gf_index_t external_A;
  gf_index_t external_B;
  double cutoff;

 public:
  // Specify return type of function to easily check type compability
  using result_t = singleton_binner;

  /// Returns integrand for the specified times
  [[nodiscard]] std::pair<result_t, int> operator()(std::vector<double> const &times) const;

  integrand_g_direct_time(g0_keldysh_contour_t g0_, gf_index_t external_A_, gf_index_t external_B_, double cutoff_ = 0.)
     : g0(std::move(g0_)), external_A(std::move(external_A_)), external_B(std::move(external_B_)), cutoff(cutoff_) {

    if ((external_A.contour.time < 0.) || (external_B.contour.time < 0.)) {
      TRIQS_RUNTIME_ERROR << "An external point has negative time.";
    }
    double const t_max = g0.get_time_max();
    if ((external_A.contour.time > t_max) || (external_B.contour.time > t_max)) {
      TRIQS_RUNTIME_ERROR << "An external point is out of g0 time window.";
    }
    if ((external_A.orbital != 0 || external_B.orbital != 0) && g0.model.make_dot_lead == false) {
      TRIQS_RUNTIME_ERROR << "Need g0_model with make_dot_lead == true to calculate off_diagonal gf.";
    }
  };
};

} // namespace keldy::impurity_oneband
