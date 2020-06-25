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

class warper_product_1d_t {
 private:
  std::vector<std::variant<warper_product_1d_simple_nointerp_t, warper_product_1d_simple_t>> warpers_dims{};

 public:
  warper_product_1d_t() = default;

  int size() { return warpers_dims.size(); }

  void emplace_back(warper_product_1d_simple_nointerp_t&& w) { warpers_dims.emplace_back(w); }
  void emplace_back(warper_product_1d_simple_t&& w) { warpers_dims.emplace_back(w); }

  // ***

  [[nodiscard]] std::pair<std::vector<double>, double> map_reverse(std::vector<double> const &li_vec) const {
    auto xi_vec = li_vec;
    double jacobian_r = 1.0;
    for (auto [i, xi] : itertools::enumerate(xi_vec)) {
      // f1 defined in ui, so evaluate jacobian_f after map
      xi = std::visit([li = xi](auto &&arg) { return arg.ui_from_li(std::vector<double>{li}); }, warpers_dims[i]).at(0);
      jacobian_r *=
         std::visit([li = xi](auto &&arg) { return arg.jacobian_reverse(std::vector<double>{li}); }, warpers_dims[i]);
    }
    return std::make_pair(xi_vec, jacobian_r);
  }

  [[nodiscard]] std::pair<std::vector<double>, double> map_forward(std::vector<double> const &ui_vec) const {
    auto xi_vec = ui_vec;
    double jacobian_f = 1.0;
    for (auto [i, xi] : itertools::enumerate(xi_vec)) {
      // f1 defined in ui, so evaluate jacobian_f before map
      jacobian_f *=
         std::visit([ui = xi](auto &&arg) { return arg.jacobian_forward(std::vector<double>{ui}); }, warpers_dims[i]);
      xi = std::visit([ui = xi](auto &&arg) { return arg.li_from_ui(std::vector<double>{ui}); }, warpers_dims[i]).at(0);
    }
    return std::make_pair(xi_vec, jacobian_f);
  }

  [[nodiscard]] std::vector<double> ui_from_li(std::vector<double> const &li_vec) const {
    std::vector<double> result = li_vec;
    for (auto [i, li] : itertools::enumerate(result)) {
      li = std::visit([li = li](auto &&arg) { return arg.ui_from_li(std::vector<double>{li}); }, warpers_dims[i]).at(0);
    }
    return result;
  }

  [[nodiscard]] std::vector<double> li_from_ui(std::vector<double> const &ui_vec) const {
    auto result = ui_vec;
    for (auto [i, ui] : itertools::enumerate(result)) {
      ui = std::visit([ui = ui](auto &&arg) { return arg.li_from_ui(std::vector<double>{ui}); }, warpers_dims[i]).at(0);
    }
    return result;
  }

  [[nodiscard]] double jacobian_reverse(std::vector<double> const &li_vec) const {
    double result = 1.0;
    for (auto const [i, li] : itertools::enumerate(li_vec)) {
      result *=
         std::visit([li = li](auto &&arg) { return arg.jacobian_reverse(std::vector<double>{li}); }, warpers_dims[i]);
    }
    return result;
  }

  [[nodiscard]] double jacobian_forward(std::vector<double> const &ui_vec) const {
    double result = 1.0;
    for (auto const [i, ui] : itertools::enumerate(ui_vec)) {
      result *=
         std::visit([ui = ui](auto &&arg) { return arg.jacobian_forward(std::vector<double>{ui}); }, warpers_dims[i]);
    }
    return result;
  }

  [[nodiscard]] double operator()(std::vector<double> const &ui_vec) const {
    double result = 1.0;
    for (auto const [i, ui] : itertools::enumerate(ui_vec)) {
      result *= std::visit([ui = ui](auto &&arg) { return arg(std::vector<double>{ui}); }, warpers_dims[i]);
    }
    return result;
  }
};

// // *****
// // Maker functions

inline warper_product_1d_t make_product_1d_inverse_cube_alternate(int order, double time, double warper_scale,
                                                                  int nr_sample_points_warper) {
  warper_product_1d_t result{};

  auto f1 = [warper_scale](double t) -> double { return warper_scale / (warper_scale + t); };
  auto f2 = [warper_scale](double t) -> double {
    return warper_scale * warper_scale * warper_scale / ((warper_scale + t) * (warper_scale + t) * (warper_scale + t));
  };
  std::vector<std::function<double(double)>> f_list = {};
  for (int n = 1; n <= order; ++n) {
    if (n % 2 == 0) {
      result.emplace_back(warper_product_1d_simple_t{f2, time, nr_sample_points_warper});
    } else {
      result.emplace_back(warper_product_1d_simple_t{f1, time, nr_sample_points_warper});
    }
  }
  return result;
}

} // namespace keldy::warpers
