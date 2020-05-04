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

namespace keldy::binner {

using namespace triqs::arrays;

class outofbinner_exception_t : public std::exception {
} out_of_binner_ex;

using double_or_int = std::variant<double, unsigned int>;
using cont_coord_t = std::variant_alternative_t<0, double_or_int>;
using disc_coord_t = std::variant_alternative_t<1, double_or_int>;

template <size_t axis, size_t N, typename Coord>
auto cast_coord(Coord var) {
  if constexpr (axis < N) {
    return static_cast<cont_coord_t>(var);
  }
  if constexpr (axis >= N) {
    return static_cast<disc_coord_t>(var);
  }
};

template <unsigned int N, unsigned int M = 0>
class sparse_binner_t {
 private:
  using coord_arr_t = std::array<double_or_int, N + M>;

  coord_arr_t tmp_coord;

  template <size_t... axes, typename... Coord>
  sparse_binner_t &operator_par_impl(std::index_sequence<axes...>, Coord... coords) {
    tmp_coord = {cast_coord<axes, N>(coords)...};
    return *this;
  };

 public:
  std::vector<std::pair<coord_arr_t, dcomplex>> data;

  template <typename T>
  auto &operator*=(T scalar) {
    for (auto &rh : data) {
      rh.second *= scalar;
    }
    return *this;
  };

  template <typename... Coord>
  sparse_binner_t &operator()(Coord... coords) {
    return operator_par_impl(std::make_index_sequence<N + M>(), coords...);
  };

  void operator<<(dcomplex value) {
    auto loc =
       std::find_if(std::begin(data), std::end(data), [coord = tmp_coord](auto &el) { return (el.first == coord); });
    if (loc != std::end(data)) {
      (*loc).second += value;
    } else {
      data.push_back(make_pair(tmp_coord, value));
    }
  };
};

template <unsigned N, unsigned M = 0>
class binner_t {
  using coord_arr_t = std::array<double_or_int, N + M>;

 private:
  array<dcomplex, N + M> data;
  std::array<std::tuple<double, double, size_t, double>, N> continuous_axes; // (xmin, xmax, nr_bins, bin_size)
  std::array<size_t, M> discreet_axes;

  array<unsigned long, N + M> nr_values_added;
  unsigned long nr_values_dropped = 0;

  coord_arr_t tmp_coord;

  template <size_t axis>
  [[nodiscard]] inline size_t coord2index(double_or_int x_) const {
    if constexpr (axis < N) {
      auto x = std::get<cont_coord_t>(x_);
      auto [xmin, xmax, n, bin_size] = continuous_axes[axis];

      if (xmin <= x && x < xmax) {
        return int((x - xmin) / bin_size);
      }
      if (x == xmax) {
        return n - 1;
      }
    } else if constexpr (axis >= N) {
      auto x = std::get<disc_coord_t>(x_);
      if (0 <= x && x < discreet_axes[axis - N]) {
        return x;
      }
    }
    throw out_of_binner_ex;
  };

  template <typename... Indices>
  inline void accumulate_impl_core(dcomplex value, Indices... indices) {
    static_assert(sizeof...(Indices) == N + M);
    data(indices...) += value;
    nr_values_added(indices...)++;
  };

  template <size_t... axes, typename... Coord>
  binner_t &operator_par_impl(std::index_sequence<axes...>, Coord... coords) {
    tmp_coord = {cast_coord<axes, N>(coords)...};
    return *this;
  };

  template <size_t... axes>
  inline void accumulate_impl(dcomplex value, std::index_sequence<axes...>, coord_arr_t coords) {
    try {
      accumulate_impl_core(value, coord2index<axes>(coords[axes])...);
    } catch (outofbinner_exception_t &e) {
      nr_values_dropped++;
    }
  };

  //template<size_t... axes, typename Coord>
  //inline void accumulate_impl(dcomplex value, std::index_sequence<axes...>, std::array<Coord, N+M> coords_){
  //  coord_arr_t coords {cast_coord<axes, N>(coords_[axes])...};
  //  accumulate_impl(value, std::make_index_sequence<N+M>(), coords);
  //};

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

  void accumulate(coord_arr_t coords, dcomplex value) {
    accumulate_impl(value, std::make_index_sequence<N + M>(), coords);
  };

  template <typename... Coord>
  binner_t &operator()(Coord... coords) {
    return operator_par_impl(std::make_index_sequence<N + M>(), coords...);
  };

  void operator<<(dcomplex value) { accumulate(tmp_coord, value); };

  [[nodiscard]] auto get_data() const { return data; };
  [[nodiscard]] auto get_nr_values_added() const { return nr_values_added; };
  [[nodiscard]] auto get_nr_values_dropped() const { return nr_values_dropped; };

  [[nodiscard]] auto get_bin_coord(size_t axis = 0) const {
    auto [xmin, xmax, n, bin_size] = continuous_axes[axis];
    array<double, 1> bin_coord(n);
    for (int i = 0; i < n; ++i) {
      bin_coord(i) = xmin + (i + 0.5) * bin_size;
    }
    return bin_coord;
  };

  [[nodiscard]] auto get_bin_size(size_t axis = 0) const { return std::get<3>(continuous_axes[axis]); };

  auto &operator+=(sparse_binner_t<N, M> const &rhs) {
    for (const auto &rh : rhs.data) {
      accumulate(rh.first, rh.second);
    }
    return *this;
  }

  template <typename T>
  auto &operator*=(T scalar) {
    data *= scalar;
    return *this;
  }
};

} // namespace keldy::binner
