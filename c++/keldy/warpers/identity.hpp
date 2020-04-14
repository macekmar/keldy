/*******************************************************************************
 *
 * keldy
 *
 * Copyright (c) 2020 The Simons Foundation
 *   authors: Philipp Dumitrescu
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
#include <vector>

namespace keldy::warpers {

class warper_identity_t {
 public:
  warper_identity_t() = default;

  [[nodiscard]] std::pair<std::vector<double>, double> map_reverse(std::vector<double> const &li_vec) const {
    return std::make_pair(ui_from_li(li_vec), jacobian_reverse(li_vec));
  }

  [[nodiscard]] std::pair<std::vector<double>, double> map_forward(std::vector<double> const &ui_vec) const {
    return std::make_pair(li_from_ui(ui_vec), jacobian_forward(ui_vec));
  }

  [[nodiscard]] std::vector<double> ui_from_li(std::vector<double> const &li_vec) const { return li_vec; }

  [[nodiscard]] std::vector<double> li_from_ui(std::vector<double> const &ui_vec) const { return ui_vec; }

  [[nodiscard]] double jacobian_reverse([[maybe_unused]] std::vector<double> const &li_vec) const { return 1.0; }

  [[nodiscard]] double jacobian_forward([[maybe_unused]] std::vector<double> const &ui_vec) const { return 1.0; }
};

} // namespace keldy::warpers
