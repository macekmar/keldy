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

using namespace triqs::arrays;

class outofbinner_exception : public std::exception {};

using cont_coord_t = double;
using disc_coord_t = long;
using double_or_int = std::variant<cont_coord_t, disc_coord_t>;

template <size_t axis, size_t N, typename Coord>
inline auto cast_coord(Coord var) {
  if constexpr (axis < N) {
    return static_cast<cont_coord_t>(var);
  }
  if constexpr (axis >= N) {
    return static_cast<disc_coord_t>(var);
  }
};

template <unsigned int N, unsigned int M = 0>
class sparse_binner_t {
  static_assert(N > 0);

 private:
  using coord_arr_t = std::array<double_or_int, N + M>;

  template <size_t... axes, typename... Coord>
  void accumulate_impl(dcomplex value, std::index_sequence<axes...>, Coord... coords) {
    coord_arr_t coord_array = {cast_coord<axes, N>(coords)...};
    auto loc =
       std::find_if(std::begin(data), std::end(data), [coord = coord_array](auto &el) { return (el.first == coord); });
    if (loc != std::end(data)) {
      (*loc).second += value;
    } else {
      data.push_back(make_pair(coord_array, value));
    }
  };

 public:
  std::vector<std::pair<coord_arr_t, dcomplex>> data;

  //template <typename T>
  auto &operator*=(dcomplex scalar) {
    for (auto &rh : data) {
      rh.second *= scalar;
    }
    return *this;
  };

  template <typename... Coord>
  void accumulate(dcomplex value, Coord... coords) {
    static_assert(sizeof...(coords) == N + M);
    accumulate_impl(value, std::make_index_sequence<N + M>(), coords...);
  };

  [[nodiscard]] friend double sum_moduli(sparse_binner_t const &in) {
    double res = 0;
    for (const auto &p : in.data) {
      res += std::abs(p.second);
    }
    return res;
  };
};

///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
template <unsigned N, unsigned M = 0>
class binner_t {
  static_assert(N > 0);
  using coord_arr_t = std::array<double_or_int, N + M>;

 private:
  array<dcomplex, N + M> data;
  std::array<std::tuple<double, double, size_t, double>, N> continuous_axes; // (xmin, xmax, nr_bins, bin_size)
  std::array<size_t, M> discreet_axes;

  array<unsigned long, N + M> nr_values_added;
  unsigned long nr_values_dropped = 0;

  template <size_t axis, typename T>
  [[nodiscard]] inline size_t coord2index(T x) const {
    if constexpr (axis < N) {
      auto [xmin, xmax, n, bin_size] = continuous_axes[axis];

      if (xmin <= x && x < xmax) {
        return int((x - xmin) / bin_size);
      }
      if (x == xmax) {
        return n - 1;
      }
      throw outofbinner_exception();
    } else if constexpr (axis >= N) {
      if (0 <= x && x < discreet_axes[axis - N]) {
        return x;
      }
      throw outofbinner_exception();
    }
  };

  template <size_t axis>
  [[nodiscard]] inline size_t coord2index(double_or_int x) const {
    if constexpr (axis < N) {
      return coord2index<axis>(std::get<cont_coord_t>(x));
    } else if constexpr (axis >= N) {
      return coord2index<axis>(std::get<disc_coord_t>(x));
    }
  };

  template <typename... Indices>
  inline void accumulate_impl_core(dcomplex value, Indices... indices) {
    static_assert(sizeof...(Indices) == N + M);
    data(indices...) += value;
    nr_values_added(indices...)++;
  };

  template <size_t... axes, typename... Coord>
  inline void accumulate_impl(dcomplex value, std::index_sequence<axes...>, Coord... coords) {
    try {
      accumulate_impl_core(value, coord2index<axes>(coords)...);
    } catch (outofbinner_exception &e) {
      nr_values_dropped++;
    }
  };

  template <size_t... axes>
  inline void accumulate_impl_array(dcomplex value, std::index_sequence<axes...>, coord_arr_t coords) {
    accumulate_impl(value, std::make_index_sequence<N + M>(), coords[axes]...);
  };

 public:
  binner_t(std::array<std::tuple<double, double, size_t>, N> _continuous_axes,
           std::array<size_t, M> _discreet_axes = {})
     : discreet_axes(std::move(_discreet_axes)) {

    triqs::utility::mini_vector<size_t, N + M> shape;
    for (auto [k, tuple] : itertools::enumerate(_continuous_axes)) {
      auto [xmin, xmax, n] = tuple;
      if (xmin >= xmax) {
        TRIQS_RUNTIME_ERROR << "Boundaries of binning range are in wrong order (" << xmin << " >= " << xmax << ")";
      }
      continuous_axes[k] = {xmin, xmax, n, (xmax - xmin) / n};
      shape[k] = n;
    }
    for (auto [k, n] : itertools::enumerate(discreet_axes)) {
      shape[N + k] = n;
    }
    data = array<dcomplex, N + M>(shape);
    data() = 0;
    nr_values_added = array<unsigned long, N + M>(shape);
    nr_values_added() = 0;
  };

  template <typename... Coord>
  void accumulate(dcomplex value, Coord... coords) {
    static_assert(sizeof...(coords) == N + M);
    accumulate_impl(value, std::make_index_sequence<N + M>(), coords...);
  };

  [[nodiscard]] auto get_data() const { return data; };
  [[nodiscard]] auto get_nr_values_added() const { return nr_values_added; };
  [[nodiscard]] auto get_nr_values_dropped() const { return nr_values_dropped; };
  [[nodiscard]] auto get_continuous_axes() const { return continuous_axes; };
  [[nodiscard]] auto get_discreet_axes() const { return discreet_axes; };

  [[nodiscard]] auto get_bin_coord(size_t axis = 0) const {
    if (axis >= N) {
      TRIQS_RUNTIME_ERROR << "Axis " << axis << " is discreet.";
    }
    if (axis >= N + M) {
      TRIQS_RUNTIME_ERROR << "Axis " << axis << " is larger than the number of dimensions.";
    }
    auto [xmin, xmax, n, bin_size] = continuous_axes[axis];
    array<double, 1> bin_coord(n);
    for (int i = 0; i < n; ++i) {
      bin_coord(i) = xmin + (i + 0.5) * bin_size;
    }
    return bin_coord;
  };

  [[nodiscard]] auto get_bin_size(size_t axis = 0) const {
    if (axis >= N) {
      TRIQS_RUNTIME_ERROR << "Axis " << axis << " is discreet.";
    }
    if (axis >= N + M) {
      TRIQS_RUNTIME_ERROR << "Axis " << axis << " is larger than the number of dimensions.";
    }

    return std::get<3>(continuous_axes[axis]);
  };

  auto &operator+=(sparse_binner_t<N, M> const &rhs) {
    for (const auto &rh : rhs.data) {
      accumulate_impl_array(rh.second, std::make_index_sequence<N + M>(), rh.first);
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

    bool all_eq = true;
    for (int ax = 0; ax < N; ++ax) {
      auto const xmin_vec = mpi::gather(std::vector<double>({std::get<0>(in.continuous_axes[ax])}), c, root, all);
      auto const xmax_vec = mpi::gather(std::vector<double>({std::get<1>(in.continuous_axes[ax])}), c, root, all);
      auto const nr_bins_vec = mpi::gather(std::vector<size_t>({std::get<2>(in.continuous_axes[ax])}), c, root, all);
      all_eq = all_eq && std::equal(xmin_vec.cbegin() + 1, xmin_vec.cend(), xmin_vec.cbegin());
      all_eq = all_eq && std::equal(xmax_vec.cbegin() + 1, xmax_vec.cend(), xmax_vec.cbegin());
      all_eq = all_eq && std::equal(nr_bins_vec.cbegin() + 1, nr_bins_vec.cend(), nr_bins_vec.cbegin());
    }
    for (int ax = N; ax < M; ++ax) {
      auto const nr_bins_vec = mpi::gather(std::vector<size_t>({in.discreet_axes[ax]}), c, root, all);
      all_eq = all_eq && std::equal(nr_bins_vec.cbegin() + 1, nr_bins_vec.cend(), nr_bins_vec.cbegin());
    }

    if (!all_eq) {
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
