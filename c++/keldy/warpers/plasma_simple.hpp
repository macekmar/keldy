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
#include "../gsl_interp_wrap.hpp"
#include <algorithm>
#include <any>
#include <functional>
#include <numeric>
#include <triqs/gfs.hpp>

// Plasma 1d -- vector of function

namespace keldy {

using namespace triqs::gfs;

// Coordinate Transforms (ui <-> vi):
inline std::vector<double> vi_from_ui(double t_max, std::vector<double> const &u_times) {
  std::vector<double> v_times = u_times;
  std::sort(v_times.begin(), v_times.end(), std::greater<>());
  int n = v_times.size();
  for (int i = n - 1; i > 0; i--) {
    v_times[i] = (v_times[i - 1] - v_times[i]);
  }
  v_times[0] = t_max - v_times[0];
  return v_times;
}

inline std::vector<double> ui_from_vi(double t_max, std::vector<double> const &v_times) {
  std::vector<double> u_times = v_times;
  u_times[0] = t_max - u_times[0];
  std::partial_sum(u_times.begin(), u_times.end(), u_times.begin(), std::minus<>());
  return u_times;
}

// Warpers:

// Todo / Qs:
// * Should we use gfs to store data. Or other?
// * replace nr_function_sample_points by variable sample_grid?
// * OP: wrapping of constructor with funciton for cpp2py. No templating?

struct idenity_function {
  double operator()(double t) { return 1.0; }
};

// std::function<dcomplex(double)>;
using gf_t = triqs::gfs::gf<retime, scalar_real_valued>;

class warper_plasma_simple_t {
 private:
  double t_max{};
  double f1_integrate_norm{};
  gf_t f1_integrated;
  gf_t f1_integrated_inverse;

  std::function<double(double)> f1;

 public:
  // might need if want pass a lambda: but should not need.
  // template <typename F>// CPP2PY Should ignore
  // warper_plasma_simple_t(F f1_, double t_max_, int nr_function_sample_points) :
  //     warper_plasma_simple_t(std::function<double(double)>(f1_), t_max_, nr_function_sample_points) {} // points vs resampling points

  // Identity Constructor: should use this if nothing else is specified
  warper_plasma_simple_t(double t_max_) : warper_plasma_simple_t{idenity_function{}, t_max_, 4} {}

  warper_plasma_simple_t(std::function<double(double)> f1_, double t_max_,
                         int nr_function_sample_points) // points vs resampling points
     : t_max(t_max_), f1(std::move(f1_)) {
    // Integrate Ansatz using Trapezoid Rule
    f1_integrated = gf_t({0.0, t_max, nr_function_sample_points});
    double delta = f1_integrated.mesh().delta();
    for (auto const &t : f1_integrated.mesh()) {
      f1_integrated[t] = f1(t - delta / 2) * delta;
    }

    auto &data = f1_integrated.data();
    data(0) = 0.0;

    std::partial_sum(data.begin(), data.end(), data.begin());
    f1_integrate_norm = data(data.size() - 1);
    data /= f1_integrate_norm; // normalize

    // Inverse Function via interpolation
    triqs::arrays::array<double, 1> mesh_time(nr_function_sample_points);
    for (auto const &[i, t] : itertools::enumerate(f1_integrated.mesh())) {
      mesh_time(i) = t;
    }

    f1_integrated_inverse = gf_t({0.0, 1.0, 1 + 5 * (nr_function_sample_points - 1)});
    details::gsl_interp_wrapper_t interpolate(gsl_interp_steffen, data, mesh_time);
    for (auto l : f1_integrated_inverse.mesh()) {
      f1_integrated_inverse[l] = interpolate(l);
    }
  }

  std::vector<double> ui_from_li(std::vector<double> const &li_vec) const {
    std::vector<double> result = li_vec;
    for (auto &li : result) {
      li = f1_integrated_inverse(li);
    }
    return ui_from_vi(t_max, result);
  }

  std::vector<double> li_from_ui(std::vector<double> const &ui_vec) const {
    auto result = vi_from_ui(t_max, ui_vec);
    // map ui to vi
    for (auto &ui : result) {
      ui = f1_integrated(ui);
    }
    return result;
  }

  double jacobian(std::vector<double> const &li_vec) const {
    double result = 1.0;
    for (auto li : li_vec) {
      result *= f1_integrate_norm / f1(f1_integrated_inverse(li));
    }
    return result;
  }

  double evaluate_warping_function(std::vector<double> const &ui_vec) const {
    double result = 1.0;
    auto vi_vec = vi_from_ui(t_max, ui_vec);
    for (auto vi : vi_vec) {
      result *= f1(vi);
    }
    return result;
  }
};

} // namespace keldy
