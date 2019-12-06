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

#include "../common.hpp"
// #include "plasma.hpp"
#include "../interfaces/gsl_interp_wrap.hpp"
#include <algorithm>
#include <any>
#include <functional>
#include <numeric>
#include <triqs/gfs.hpp>

// Plasma 1d -- vector of functions

namespace keldy {

using namespace triqs::gfs;

// Warpers:

// Todo / Qs:
// * Should we use gfs to store data. Or other?
// * replace nr_function_sample_points by variable sample_grid?
// * OP: wrapping of constructor with funciton for cpp2py. No templating?

// std::function<dcomplex(double)>;
using gf_t = triqs::gfs::gf<retime, scalar_real_valued>;

class warper_plasma_1D_t {
 private:
  double t_max{};
  std::vector<double> fn_integrate_norm{};
  std::vector<gf_t> fn_integrated;
  std::vector<gf_t> fn_integrated_inverse;


  std::vector<std::function<double(double)>> fn;

 public:
  
    // Identity Constructor: should use this if nothing else is specified
  warper_plasma_1D_t(double t_max_) : warper_plasma_1D_t{{idenity_function{}}, t_max_, 8} {}

  warper_plasma_1D_t(std::vector<std::function<double(double)>> fn_, double t_max_,
                         int nr_function_sample_points) // points vs resampling points
     : t_max(t_max_), fn(std::move(fn_)) {
    for (int axis = 0; axis < fn.size(); axis++) {
      // Integrate Ansatz using Trapezoid Rule
      fn_integrated.push_back(gf_t({0.0, t_max, nr_function_sample_points}));
      double delta = fn_integrated[axis].mesh().delta();
      for (auto const &t : fn_integrated[axis].mesh()) {
        fn_integrated[axis][t] = fn[axis](t - delta / 2) * delta;
      }
      auto &data = fn_integrated[axis].data();
      data(0) = 0.0;

      std::partial_sum(data.begin(), data.end(), data.begin());
      fn_integrate_norm.push_back(data(data.size() - 1));
      data /= fn_integrate_norm[axis]; // normalize each axis separtely
    
    
      // Inverse Function via interpolation
      triqs::arrays::array<double, 1> mesh_time(nr_function_sample_points);
      for (auto const &[i, t] : itertools::enumerate(fn_integrated[axis].mesh())) {
        mesh_time(i) = t;
      }

      fn_integrated_inverse.push_back(gf_t({0.0, 1.0, 1 + 5 * (nr_function_sample_points - 1)}));
      details::gsl_interp_wrapper_t interpolate(gsl_interp_akima, data, mesh_time); //gsl_interp_steffen
      for (auto l : fn_integrated_inverse[axis].mesh()) {
        fn_integrated_inverse[axis][l] = interpolate(l);
      }
    }
  }

  std::vector<double> ui_from_li(std::vector<double> const &li_vec) const {
    std::vector<double> result = li_vec;
    for (auto [i, li] : itertools::enumerate(result)) {
      li = fn_integrated_inverse[i](li);
    }
    return ui_from_vi(t_max, result);
  }

  std::vector<double> li_from_ui(std::vector<double> const &ui_vec) const {
    auto result = vi_from_ui(t_max, ui_vec);
    // map ui to vi
    for (auto [i, ui] : itertools::enumerate(result)) {
      ui = fn_integrated[i](ui);
    }
    return result;
  }

  double jacobian(std::vector<double> const &li_vec) const {
    double result = 1.0;
    for (auto [i, li] : itertools::enumerate(li_vec)) {
      result *= fn_integrate_norm[i] / fn[i](fn_integrated_inverse[i](li));
    }
    return result;
  }

  double evaluate_warping_function(std::vector<double> const &ui_vec) const {
    double result = 1.0;
    auto vi_vec = vi_from_ui(t_max, ui_vec);
    for (auto [i,vi] : itertools::enumerate(vi_vec)) {
      result *= fn[i](vi);
    }
    return result;
  }
};

} // namespace keldy
