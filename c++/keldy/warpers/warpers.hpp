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

// // warper eraser
// class warper_t {
//  private:
//   std::any warp;
//   std::function<std::vector<double>(std::vector<double> const &)> ui_from_li_;
//   std::function<double(std::vector<double> const &)> jacobian_;

//  public:

//   warper_t() = default;

//   template <typename T>
//   CPP2PY_IGNORE warper_t(T x) : warp{std::move(x)} {
//     T const *p = std::any_cast<T const>(&warp); // retrieve a pointer on the stored T object

//     ui_from_li_ = [p](std::vector<double> const &li_vec) { return p->ui_from_li(li_vec); };
//     jacobian_  = [p](std::vector<double> const &li_vec) { return p->jacobian(li_vec); };
//   }

//   double jacobian(std::vector<double> const &x) const {return jacobian_(x);}
//   std::vector<double> ui_from_li(std::vector<double> const &x) const {return ui_from_li_(x);}

// };

} // namespace keldy
