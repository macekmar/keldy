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

#include "warpers_common.hpp"
#include "identity.hpp"
#include "plasma_uv.hpp"
#include "product_1d_simple.hpp"
#include "product_1d.hpp"
#include "plasma_projection.hpp"

#include <variant>
#include <vector>

namespace keldy {

// varient is default constructable to hold value of first alternative (if that is default constructable)
using warper_variant =
   std::variant<warper_identity_t, warper_plasma_uv_t, warper_product_1d_simple_t, warper_product_1d_t, warper_plasma_projection_t>;

class warper_train_t {
 public:
  std::vector<warper_variant> warpers{};

  std::vector<double> ui_from_li(std::vector<double> const &li_vec) const {
    std::vector<double> result = li_vec;
    for (auto it = warpers.crbegin(); it != warpers.crend(); it++) { // REVERSE
      result = std::visit([&result](auto &&arg) -> std::vector<double> { return arg.ui_from_li(result); }, *it);
    }
    return result;
  }

  std::vector<double> li_from_ui(std::vector<double> const &ui_vec) const {
    auto result = ui_vec;
    for (auto &v : warpers) { // FORWARD
      result = std::visit([&result](auto &&arg) -> std::vector<double> { return arg.li_from_ui(result); }, v);
    }
    return result;
  }

  double jacobian(std::vector<double> const &li_vec) const {
    double result = 1.0;
    std::vector<double> li_vec_tmp = li_vec;
    for (auto it = warpers.crbegin(); it != warpers.crend(); it++) { // REVERSE
      result *= std::visit([&li_vec_tmp](auto &&arg) { return arg.jacobian(li_vec_tmp); }, *it);
      li_vec_tmp = std::visit([&li_vec_tmp](auto &&arg) { return arg.ui_from_li(li_vec_tmp); }, *it);
    }
    return result;
  }

  double operator()(std::vector<double> const &ui_vec) const {
    double result = 1.0;
    std::vector<double> ui_vec_tmp = ui_vec;
    for (auto &v : warpers) { // FORWARD
      result *= std::visit([&ui_vec_tmp](auto &&arg) -> double { return arg(ui_vec_tmp); }, v);
      ui_vec_tmp = std::visit([&ui_vec_tmp](auto &&arg) { return arg.li_from_ui(ui_vec_tmp); }, v);
    }
    return result;
  }
};

} // namespace keldy