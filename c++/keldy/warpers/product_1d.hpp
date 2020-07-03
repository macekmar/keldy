/*******************************************************************************
 *
 * keldy
 *
 * Copyright (c) 2020 The Simons Foundation
 * Copyright (c) 2020 CEA: Commissariat à l’énergie atomique
 *                         et aux énergies alternatives
 *   authors: Philipp Dumitrescu, Marjan Macek
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

#include "../common.hpp"
#include "warpers_common.hpp"
#include "product_1d_simple.hpp"
#include "../interfaces/gsl_interp_wrap.hpp"
#include <algorithm>
#include <any>
#include <functional>
#include <numeric>
#include <triqs/gfs.hpp>

namespace keldy::warpers {

using gf_t = triqs::gfs::gf<triqs::gfs::retime, triqs::gfs::scalar_real_valued>;

// r standalone functions, as to not template classes (for python wrapping)

template <typename W>
[[nodiscard]] std::pair<std::vector<double>, double> warper_1d_map_reverse(std::vector<W> const &warpers_dims,
                                                                           std::vector<double> const &li_vec) {
  auto xi_vec = li_vec;
  double jacobian_r = 1.0;

  for (int i = 0; i < xi_vec.size(); i++) {
    auto &xi = xi_vec[i];
    auto [xi_vec_1, jac_r_1] = warpers_dims[i].map_reverse(std::vector<double>{xi});
    xi = xi_vec_1.at(0);
    jacobian_r *= jac_r_1;
  }
  return std::make_pair(xi_vec, jacobian_r);
}

template <typename W>
[[nodiscard]] std::pair<std::vector<double>, double> warper_1d_map_forward(std::vector<W> const &warpers_dims,
                                                                           std::vector<double> const &ui_vec) {
  auto xi_vec = ui_vec;
  double jacobian_f = 1.0;
  for (int i = 0; i < xi_vec.size(); i++) {
    auto &xi = xi_vec[i];
    auto [xi_vec_1, jac_f_1] = warpers_dims[i].map_forward(std::vector<double>{xi});
    xi = xi_vec_1.at(0);
    jacobian_f *= jac_f_1;
  }
  return std::make_pair(xi_vec, jacobian_f);
}

template <typename W>
[[nodiscard]] std::vector<double> warper_1d_ui_from_li(std::vector<W> const &warpers_dims,
                                                       std::vector<double> const &li_vec) {
  std::vector<double> result = li_vec;
  for (int i = 0; i < result.size(); i++) {
    auto &li = result[i];
    li = warpers_dims[i].ui_from_li(std::vector<double>{li}).at(0);
  }
  return result;
}

template <typename W>
[[nodiscard]] std::vector<double> warper_1d_li_from_ui(std::vector<W> const &warpers_dims,
                                                       std::vector<double> const &ui_vec) {
  auto result = ui_vec;
  for (int i = 0; i < result.size(); i++) {
    auto &ui = result[i];
    ui = warpers_dims[i].li_from_ui(std::vector<double>{ui}).at(0);
  }
  return result;
}

template <typename W>
[[nodiscard]] double warper_1d_jacobian_reverse(std::vector<W> const &warpers_dims, std::vector<double> const &li_vec) {
  double result = 1.0;
  for (auto const [i, li] : itertools::enumerate(li_vec)) {
    result *= warpers_dims[i].jacobian_reverse(std::vector<double>{li});
  }
  return result;
}

template <typename W>
[[nodiscard]] double warper_1d_jacobian_forward(std::vector<W> const &warpers_dims, std::vector<double> const &ui_vec) {
  double result = 1.0;
  for (auto const [i, ui] : itertools::enumerate(ui_vec)) {
    result *= warpers_dims[i].jacobian_forward(std::vector<double>{ui});
  }
  return result;
}

template <typename W>
[[nodiscard]] double warper_1d_evalutate(std::vector<W> const &warpers_dims, std::vector<double> const &ui_vec) {
  double result = 1.0;
  for (auto const [i, ui] : itertools::enumerate(ui_vec)) {
    result *= warpers_dims[i](std::vector<double>{ui});
  }
  return result;
}

// ****************************************************************************

class warper_product_1d_t {
 private:
  std::vector<warper_product_1d_simple_t> warpers_dims{};

 public:
  warper_product_1d_t() = default;

  int size() { return warpers_dims.size(); }
  void emplace_back(warper_product_1d_simple_t w) { warpers_dims.emplace_back(std::move(w)); }

  [[nodiscard]] std::pair<std::vector<double>, double> map_reverse(std::vector<double> const &li_vec) const {
    return warper_1d_map_reverse(warpers_dims, li_vec);
  }

  [[nodiscard]] std::pair<std::vector<double>, double> map_forward(std::vector<double> const &ui_vec) const {
    return warper_1d_map_forward(warpers_dims, ui_vec);
  }

  [[nodiscard]] std::vector<double> ui_from_li(std::vector<double> const &li_vec) const {
    return warper_1d_ui_from_li(warpers_dims, li_vec);
  }

  [[nodiscard]] std::vector<double> li_from_ui(std::vector<double> const &ui_vec) const {
    return warper_1d_li_from_ui(warpers_dims, ui_vec);
  }

  [[nodiscard]] double jacobian_reverse(std::vector<double> const &li_vec) const {
    return warper_1d_jacobian_reverse(warpers_dims, li_vec);
  }

  [[nodiscard]] double jacobian_forward(std::vector<double> const &ui_vec) const {
    return warper_1d_jacobian_forward(warpers_dims, ui_vec);
  }

  [[nodiscard]] double operator()(std::vector<double> const &ui_vec) const {
    return warper_1d_evalutate(warpers_dims, ui_vec);
  }
};

class warper_product_1d_interp_nearest_t {
 private:
  std::vector<warper_product_1d_simple_interp_nearest_t> warpers_dims{};

 public:
  warper_product_1d_interp_nearest_t() = default;

  int size() { return warpers_dims.size(); }
  void emplace_back(warper_product_1d_simple_interp_nearest_t w) { warpers_dims.emplace_back(std::move(w)); }

  [[nodiscard]] std::pair<std::vector<double>, double> map_reverse(std::vector<double> const &li_vec) const {
    return warper_1d_map_reverse(warpers_dims, li_vec);
  }

  [[nodiscard]] std::pair<std::vector<double>, double> map_forward(std::vector<double> const &ui_vec) const {
    return warper_1d_map_forward(warpers_dims, ui_vec);
  }

  [[nodiscard]] std::vector<double> ui_from_li(std::vector<double> const &li_vec) const {
    return warper_1d_ui_from_li(warpers_dims, li_vec);
  }

  [[nodiscard]] std::vector<double> li_from_ui(std::vector<double> const &ui_vec) const {
    return warper_1d_li_from_ui(warpers_dims, ui_vec);
  }

  [[nodiscard]] double jacobian_reverse(std::vector<double> const &li_vec) const {
    return warper_1d_jacobian_reverse(warpers_dims, li_vec);
  }

  [[nodiscard]] double jacobian_forward(std::vector<double> const &ui_vec) const {
    return warper_1d_jacobian_forward(warpers_dims, ui_vec);
  }

  [[nodiscard]] double operator()(std::vector<double> const &ui_vec) const {
    return warper_1d_evalutate(warpers_dims, ui_vec);
  }
};

class warper_product_1d_interp_hybrid_t {
 private:
  std::vector<warper_product_1d_simple_interp_hybrid_t> warpers_dims{};

 public:
  warper_product_1d_interp_hybrid_t() = default;

  int size() { return warpers_dims.size(); }
  void emplace_back(warper_product_1d_simple_interp_hybrid_t w) { warpers_dims.emplace_back(std::move(w)); }

  [[nodiscard]] std::pair<std::vector<double>, double> map_reverse(std::vector<double> const &li_vec) const {
    return warper_1d_map_reverse(warpers_dims, li_vec);
  }

  [[nodiscard]] std::pair<std::vector<double>, double> map_forward(std::vector<double> const &ui_vec) const {
    return warper_1d_map_forward(warpers_dims, ui_vec);
  }

  [[nodiscard]] std::vector<double> ui_from_li(std::vector<double> const &li_vec) const {
    return warper_1d_ui_from_li(warpers_dims, li_vec);
  }

  [[nodiscard]] std::vector<double> li_from_ui(std::vector<double> const &ui_vec) const {
    return warper_1d_li_from_ui(warpers_dims, ui_vec);
  }

  [[nodiscard]] double jacobian_reverse(std::vector<double> const &li_vec) const {
    return warper_1d_jacobian_reverse(warpers_dims, li_vec);
  }

  [[nodiscard]] double jacobian_forward(std::vector<double> const &ui_vec) const {
    return warper_1d_jacobian_forward(warpers_dims, ui_vec);
  }

  [[nodiscard]] double operator()(std::vector<double> const &ui_vec) const {
    return warper_1d_evalutate(warpers_dims, ui_vec);
  }
};
// ****************************************************************************
// Maker functions

inline warper_product_1d_t make_product_1d_inverse_cube_alternate(int order, double time, double warper_scale) {
  warper_product_1d_t result{};
  for (int n = 1; n <= order; ++n) {
    if (n % 2 == 0) {
      result.emplace_back(make_product_1d_simple_inverse_power(3, time, warper_scale));
    } else {
      result.emplace_back(make_product_1d_simple_inverse(time, warper_scale));
    }
  }
  return result;
}

inline warper_product_1d_interp_nearest_t
make_product_1d_inverse_cube_alternate_interp_nearest(int order, double time, double warper_scale,
                                                      int nr_sample_points_warper) {
  warper_product_1d_interp_nearest_t result{};

  auto f1 = [warper_scale](double t) -> double { return warper_scale / (warper_scale + t); };
  auto f2 = [warper_scale](double t) -> double {
    return warper_scale * warper_scale * warper_scale / ((warper_scale + t) * (warper_scale + t) * (warper_scale + t));
  };
  std::vector<std::function<double(double)>> f_list = {};
  for (int n = 1; n <= order; ++n) {
    if (n % 2 == 0) {
      result.emplace_back(warper_product_1d_simple_interp_nearest_t{f2, time, nr_sample_points_warper});
    } else {
      result.emplace_back(warper_product_1d_simple_interp_nearest_t{f1, time, nr_sample_points_warper});
    }
  }
  return result;
}

inline warper_product_1d_interp_hybrid_t
make_product_1d_inverse_cube_alternate_interp_hybrid(int order, double time, double warper_scale,
                                                     int nr_sample_points_warper) {
  warper_product_1d_interp_hybrid_t result{};

  auto f1 = [warper_scale](double t) -> double { return warper_scale / (warper_scale + t); };
  auto f2 = [warper_scale](double t) -> double {
    return warper_scale * warper_scale * warper_scale / ((warper_scale + t) * (warper_scale + t) * (warper_scale + t));
  };
  std::vector<std::function<double(double)>> f_list = {};
  for (int n = 1; n <= order; ++n) {
    if (n % 2 == 0) {
      result.emplace_back(warper_product_1d_simple_interp_hybrid_t{f2, time, nr_sample_points_warper});
    } else {
      result.emplace_back(warper_product_1d_simple_interp_hybrid_t{f1, time, nr_sample_points_warper});
    }
  }
  return result;
}

} // namespace keldy::warpers
