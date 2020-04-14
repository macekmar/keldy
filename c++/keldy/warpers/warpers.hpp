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

#include "warpers_common.hpp"
#include "identity.hpp"
#include "plasma_uv.hpp"
#include "product_1d_simple.hpp"
#include "product_1d.hpp"
#include "plasma_projection.hpp"

#include <variant>
#include <vector>

namespace keldy::warpers {

// varient is default constructable to hold value of first alternative (if that is default constructable)
using warper_variant =
   std::variant<warper_identity_t, warper_plasma_uv_t, warper_product_1d_simple_t, warper_product_1d_t>;

class warper_train_t {
 public:
  std::vector<warper_variant> warpers{};

  // Joint Evolution: Coordinate-Transform and Jacobian including partial evolution
  std::pair<std::vector<double>, double> map_reverse(std::vector<double> const &li_vec, int start_domain_nr,
                                                        int end_domain_nr) const {
    // EXPECTS(0 <= start_domain_nr <= warper.size(), 0 <= end_domain_nr <= warper.size(),
    //            start_domain_nr <= end_domain_nr  )
    double jacobian_reverse_result = 1.0;
    std::vector<double> ui_vec = li_vec;

    for (auto it = warpers.crbegin() + warpers.size() - start_domain_nr; it != warpers.crend() - end_domain_nr; it++) {
      auto out = std::visit([&ui_vec](auto &&arg) { return arg.map_reverse(ui_vec); }, *it);
      ui_vec = out.first;
      jacobian_reverse_result *= out.second;
    }
    return std::make_pair(ui_vec, jacobian_reverse_result);
  }

  // better to try itertools / at??
  std::pair<std::vector<double>, double> map_forward(std::vector<double> const &ui_vec, int start_domain_nr,
                                                        int end_domain_nr) const {
    // EXPECTS(0 <= start_domain_nr <= warper.size(), 0 <= end_domain_nr <= warper.size(),
    //            start_domain_nr <= end_domain_nr  )
    double jacobian_forward_result = 1.0;
    std::vector<double> li_vec = ui_vec;

    for (auto it = warpers.cbegin() + start_domain_nr; it != warpers.cbegin() + end_domain_nr; it++) {
      auto out = std::visit([&li_vec](auto &&arg) { return arg.map_forward(li_vec); }, *it);
      li_vec = out.first;
      jacobian_forward_result *= out.second;
    }
    return std::make_pair(ui_vec, jacobian_forward_result);
  }

  // Functions which are concept of warper:
  std::pair<std::vector<double>, double> map_reverse(std::vector<double> const &li_vec) const {
    return map_reverse(li_vec, warpers.size(), 0);
  }

  std::pair<std::vector<double>, double> map_forward(std::vector<double> const &ui_vec) const {
    return map_forward(ui_vec, 0, warpers.size());
  }

  double jacobian_reverse(std::vector<double> const &li_vec) const {
    auto [ui_vec, jacobian_r] = map_reverse(li_vec);
    return jacobian_r;
  }

  double jacobian_forward(std::vector<double> const &ui_vec) const {
    auto [li_vec, jacobian_f] = map_forward(ui_vec);
    return jacobian_f;
  }

  std::vector<double> ui_from_li(std::vector<double> const &li_vec) const {
    std::vector<double> result = li_vec;
    for (auto it = warpers.crbegin(); it != warpers.crend(); it++) { // REVERSE
      result = std::visit([&result](auto &&arg) -> std::vector<double> { return arg.ui_from_li(result); }, *it);
    }
    return result;
  }

  std::vector<double> li_from_ui(std::vector<double> const &ui_vec) const {
    auto result = ui_vec;
    for (auto it = warpers.cbegin(); it != warpers.cend(); it++) { // FORWARD
      result = std::visit([&result](auto &&arg) -> std::vector<double> { return arg.li_from_ui(result); }, *it);
    }
    return result;
  }
  
};

} // namespace keldy::warpers