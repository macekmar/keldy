/*******************************************************************************
 *
 * keldy
 *
 * Copyright (c) 2019 The Simons Foundation
 *   authors: Philipp Dumitrescu, Corentin Bertrand
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

#pragma once

#include "common.hpp"
#include <triqs/arrays.hpp>
#include <itertools/itertools.hpp>
#include <exception>
#include <mpi/mpi.hpp>
#include <mpi/vector.hpp>
#include <triqs/arrays/mpi.hpp>

namespace keldy::binner {

namespace mda = triqs::arrays;

using cont_coord_t = double;
using disc_coord_t = long;

/// MPI broadcast for std::array
template <typename T, size_t N>
void mpi_broadcast(std::array<T, N> &v, mpi::communicator c = {}, int root = 0) {
  if constexpr (mpi::has_mpi_type<T>) {
    if (N != 0) MPI_Bcast(v.data(), N, mpi::mpi_type<T>::get(), root, c.get());
  } else {
    for (auto &x : v)
      mpi::broadcast(x, c, root);
  }
}

/*
 * Multidimentional sparse binner with N continuous coordinates and M discreet ones.
 * Continuous coordinates are before discreet ones.
 */
template <int N, int M = 0>
class sparse_binner_t {
  static_assert(N > 0);

 private:
  template <size_t... cont_axes, size_t... disc_axes, typename... Coord>
  inline void accumulate_impl(dcomplex value, std::index_sequence<cont_axes...>, std::index_sequence<disc_axes...>,
                              Coord... coords_) {
    auto coords = std::tie(coords_...);
    coord_arr_t coord_array = {{std::get<cont_axes>(coords)...}, {std::get<N + disc_axes>(coords)...}};

    data.push_back(make_pair(coord_array, value));
  };

 public:
  using coord_arr_t = std::pair<std::array<cont_coord_t, N>, std::array<disc_coord_t, M>>;

  /// list of stored values
  std::vector<std::pair<coord_arr_t, dcomplex>> data;

  // to multiply with jacobian e.g.
  auto &operator*=(dcomplex scalar) {
    for (auto &rh : data) {
      rh.second *= scalar;
    }
    return *this;
  };

  /// Append a value into the list `data`.
  template <typename... Coord>
  void accumulate(dcomplex value, Coord... coords) {
    static_assert(sizeof...(coords) == N + M);
    accumulate_impl(value, std::make_index_sequence<N>(), std::make_index_sequence<M>(), coords...);
  };
};

/// sum moduli of values stored in sparese binner.
//  Values with same coordinates are summed before taking the modulus.
template <int N, int M>
[[nodiscard]] double sum_moduli(sparse_binner_t<N, M> const &sp_bin) {

  std::vector<std::pair<typename sparse_binner_t<N, M>::coord_arr_t, dcomplex>> data_reduced;

  for (auto const &p : sp_bin.data) {
    auto loc =
       std::find_if(std::begin(data_reduced), std::end(data_reduced), [&p](auto &el) { return (el.first == p.first); });
    if (loc != std::end(data_reduced)) {
      (*loc).second += p.second;
    } else {
      data_reduced.push_back(p);
    }
  }

  double res = 0;
  for (auto const &q : data_reduced) {
    res += std::abs(q.second);
  }
  return res;
};

///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

struct continuous_axis_t {
  double xmin;
  double xmax;
  long nr_bins;
  double bin_size;
};

void mpi_broadcast(continuous_axis_t &in, mpi::communicator c = {}, int root = 0);

bool operator==(continuous_axis_t const &a, continuous_axis_t const &b);

struct discreet_axis_t {
  long nr_bins;
};

void mpi_broadcast(discreet_axis_t &in, mpi::communicator c = {}, int root = 0);

/*
 * Multidimentional binner with N continuous coordinates and M discreet ones.
 * Continuous coordinates are before discreet ones.
 */
template <int N, int M = 0>
class binner_t {
  static_assert(N > 0);
  using coord_arr_t = std::pair<std::array<cont_coord_t, N>, std::array<disc_coord_t, M>>;

 private:
  mda::array<dcomplex, N + M> data;
  mda::array_view<dcomplex, 1> data_view;
  std::array<continuous_axis_t, N> continuous_axes;
  std::array<discreet_axis_t, M> discreet_axes;

  mda::array<long, N + M> nr_values_added;
  mda::array_view<long, 1> nr_value_added_view;
  long nr_values_dropped = 0;

  /// Find bin index from a coordinate on a given axis.
  template <size_t axis, typename T>
  [[nodiscard]] long coord2index(T x) const {
    if constexpr (axis < N) {
      auto const &ax = continuous_axes[axis];
      if (x == ax.xmax) {
        return ax.nr_bins - 1;
      }
      return int((x - ax.xmin) / ax.bin_size);
    } else { // discreet axis
      return x;
    }
  };

  /// Find flattened index from complete list of coordinates
  template <size_t... axes, typename... Coord>
  [[nodiscard]] long coord2flatidx(std::index_sequence<axes...>, Coord... coords) {
    return data.indexmap()(coord2index<axes>(coords)...);
  };

  /// Check if coordinate lies inside the binner on a given axis.
  template <size_t axis, typename T>
  bool in_bounds(T x) {
    if constexpr (axis < N) {
      auto ax = continuous_axes[axis];
      return (ax.xmin <= x && x <= ax.xmax);
    }
    return (0 <= x && x < discreet_axes[axis - N].nr_bins);
  };

  /// Check if complete list of coordinates lies inside the binner.
  template <size_t... axes, typename... Coord>
  [[nodiscard]] long in_bounds(std::index_sequence<axes...>, Coord... coords) {
    auto all = [](auto... args) -> bool { return (... && args); };
    return all(in_bounds<axes>(coords)...);
  };

  /// Accumulate a value from a sparse binner.
  template <size_t... cont_axes, size_t... disc_axes>
  void accumulate_array(dcomplex value, std::index_sequence<cont_axes...>, std::index_sequence<disc_axes...>,
                        coord_arr_t coords) {
    (*this)(coords.first[cont_axes]..., coords.second[disc_axes]...) << value;
  };

 public:
  /*
   * Constructs a multidimentional binner.
   *
   * Arguments:
   *    continous_axes_: list of tuples (xmin, xmax, nr_bins) for each continuous axis
   *    discreet_axes_: list of number of bins for each discreet axis
   */
  binner_t(std::array<std::tuple<double, double, long>, N> continuous_axes_, std::array<long, M> discreet_axes_ = {}) {

    triqs::utility::mini_vector<int, N + M> shape;
    for (auto [k, tuple] : itertools::enumerate(continuous_axes_)) {
      auto [xmin, xmax, n] = tuple;
      if (xmin >= xmax) {
        TRIQS_RUNTIME_ERROR << "Boundaries of binning range are in wrong order (" << xmin << " >= " << xmax << ")";
      }
      continuous_axes[k] = {xmin, xmax, n, (xmax - xmin) / n};
      shape[k] = n;
    }
    for (auto [k, n] : itertools::enumerate(discreet_axes_)) {
      discreet_axes[k].nr_bins = n;
      shape[N + k] = n;
    }

    data = mda::array<dcomplex, N + M>(shape);
    data() = 0;
    data_view.rebind(mda::array_view<dcomplex, 1>({mda::make_shape(data.size())}, data.storage()));
    nr_values_added = mda::array<long, N + M>(shape);
    nr_values_added() = 0;
    nr_value_added_view.rebind(
       mda::array_view<long, 1>({mda::make_shape(nr_values_added.size())}, nr_values_added.storage()));
  };

  /// Small object containing ref to a binner_t and to a particular bin
  //  Used to produce nice accumulate syntax:
  //  binner(coord1, coord2, ...) << value;
  template <typename T, typename... Coord>
  struct accessor_t {
    T &binner;
    long flat_idx; // < 0 means out of binner

    /// Accumulate a value in the binner at given coordinates
    void operator<<(dcomplex value) {
      if (flat_idx < 0) {
        binner.nr_values_dropped++;
      } else {
        binner.data_view(flat_idx) += value;
        binner.nr_value_added_view(flat_idx)++;
      }
    }
  };

  /// nice syntax for accumulate:
  //  binner(coord1, coord2, ...) << value;
  template <typename... Coord>
  accessor_t<binner_t<N, M>, Coord...> operator()(Coord... coords) {
    static_assert(sizeof...(coords) == N + M);
    if (in_bounds(std::make_index_sequence<N + M>(), coords...)) {
      return {*this, coord2flatidx(std::make_index_sequence<N + M>(), coords...)};
    }
    return {*this, -1};
  };

  [[nodiscard]] auto const &get_data() const { return data; };
  [[nodiscard]] auto const &get_nr_values_added() const { return nr_values_added; };
  [[nodiscard]] auto const &get_nr_values_dropped() const { return nr_values_dropped; };
  [[nodiscard]] auto const &get_continuous_axes() const { return continuous_axes; };
  [[nodiscard]] auto const &get_discreet_axes() const { return discreet_axes; };

  /// Get the coordinates of bins along a given (continuous) axis.
  // Returns a 1D triqs array
  [[nodiscard]] auto get_bin_coord(int axis = 0) const {
    if (axis >= N) {
      TRIQS_RUNTIME_ERROR << "Axis " << axis << " is discreet.";
    }
    if (axis >= N + M) {
      TRIQS_RUNTIME_ERROR << "Axis " << axis << " is larger than the number of dimensions.";
    }
    auto ax = continuous_axes[axis];
    mda::array<double, 1> bin_coord(ax.nr_bins);
    for (int i = 0; i < ax.nr_bins; ++i) {
      bin_coord(i) = ax.xmin + (i + 0.5) * ax.bin_size;
    }
    return bin_coord;
  };

  /// Get the width of a bin on a given (continuous) axis.
  [[nodiscard]] double get_bin_size(int axis = 0) const {
    if (axis >= N) {
      TRIQS_RUNTIME_ERROR << "Axis " << axis << " is discreet.";
    }
    if (axis >= N + M) {
      TRIQS_RUNTIME_ERROR << "Axis " << axis << " is larger than the number of dimensions.";
    }

    return continuous_axes[axis].bin_size;
  };

  /// Accumulate a sparse binner into this binner
  auto &operator+=(sparse_binner_t<N, M> const &rhs) {
    for (const auto &rh : rhs.data) {
      accumulate_array(rh.second, std::make_index_sequence<N>(), std::make_index_sequence<M>(), rh.first);
    }
    return *this;
  }

  template <typename dcomplex>
  auto &operator*=(dcomplex scalar) {
    data *= scalar;
    return *this;
  }

  friend inline binner_t<N, M> CPP2PY_IGNORE mpi_reduce(binner_t<N, M> const &in, mpi::communicator c = {},
                                                        int root = 0, bool all = false, MPI_Op op = MPI_SUM) {

    if (op != MPI_SUM) {
      TRIQS_RUNTIME_ERROR << "mpi_reduce of binner_t can only be performed with op = MPI_SUM";
    }

    std::array<continuous_axis_t, N> continuous_axes_ref;
    std::array<discreet_axis_t, M> discreet_axes_ref;
    if (c.rank() == root) {
      continuous_axes_ref = in.continuous_axes;
      discreet_axes_ref = in.discreet_axes;
    }
    mpi::broadcast(continuous_axes_ref, c, root);
    if (M > 0) {
      mpi::broadcast(discreet_axes_ref, c, root);
    }

    if (in.continuous_axes != continuous_axes_ref || (M > 0 && in.continuous_axes != continuous_axes_ref)) {
      TRIQS_RUNTIME_ERROR << "Binners to mpi_reduce are not defined with same bins. This is not supported!";
    }

    binner_t out = in; // copy

    out.data() = mpi::reduce(in.data, c, root, all, MPI_SUM);
    out.nr_values_added() = mpi::reduce(in.nr_values_added, c, root, all, MPI_SUM);
    out.nr_values_dropped = mpi::reduce(in.nr_values_dropped, c, root, all, MPI_SUM);

    return out;
  };
};

} // namespace keldy::binner
