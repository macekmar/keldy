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

#include "keldy/common.hpp"
#include "plasma.hpp"
#include "plasma_simple.hpp"

#include <algorithm>
#include <functional>

namespace keldy {

// warper eraser
struct CPP2PY_IGNORE warper_t {
  std::any warp;
  std::function<std::vector<double>(std::vector<double> const &)> ui_from_li;
  std::function<double(std::vector<double> const &)> jacobian;

  template <typename T>
  warper_t(T &&x) : warp{std::forward<T>(x)} {
    T *p = std::any_cast<T>(&warp); // retrieve a pointer on the stored T object
    ui_from_li = [p](std::vector<double> const &li_vec) { p->ui_from_li(li_vec); };
    jacobian = [p](std::vector<double> const &li_vec) { p->jacobian(li_vec); };
  }
};

}