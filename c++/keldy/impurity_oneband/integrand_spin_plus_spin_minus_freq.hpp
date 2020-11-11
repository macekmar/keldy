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

#include "../common.hpp"
#include "../binner.hpp"
#include "model.hpp"

#include <mpi/mpi.hpp>
#include <mpi/vector.hpp>
#include <triqs/arrays/math_functions.hpp>
#include <triqs/arrays/mpi.hpp>
#include <triqs/utility/exceptions.hpp>
#include <triqs/utility/first_include.hpp>
#include <algorithm>
#include <vector>
#include <tuple>

namespace keldy::impurity_oneband {

using namespace triqs::arrays;
using namespace triqs::gfs;

class chi_function_t {
  mda::array<double, 1> const omegas;
  double const time;
  gf<retime, tensor_valued<3>> chi;

 public:
  chi_function_t(g0_model_omega &model_omega, mda::array<double, 1> omegas_, double time_, keldysh_idx_t a)
     : omegas(std::move(omegas_)), time(time_) {
    auto param = model_omega.get_param();
    int const N = omegas.size();

    auto time_mesh = gf_mesh<retime>({-param.time_max, param.time_max, param.nr_time_points_gf});
    auto omega_mesh = make_adjoint_mesh(time_mesh);
    gf<refreq, tensor_valued<3>> chi_omega{omega_mesh, {N, 2, 2}}; // Keldysh matrix, not supporting orbitals yet
    for (auto [i, w0] : itertools::enumerate(omegas)) {

      for (auto w : omega_mesh) {
        chi_omega[w](i, 0, 0) =
           model_omega.g0_dot_keldysh(a, forward, w + w0) * model_omega.g0_dot_keldysh(forward, a, w + 0.);
        chi_omega[w](i, 0, 1) =
           model_omega.g0_dot_keldysh(a, forward, w + w0) * model_omega.g0_dot_keldysh(backward, a, w + 0.);
        chi_omega[w](i, 1, 0) =
           model_omega.g0_dot_keldysh(a, backward, w + w0) * model_omega.g0_dot_keldysh(forward, a, w + 0.);
        chi_omega[w](i, 1, 1) =
           model_omega.g0_dot_keldysh(a, backward, w + w0) * model_omega.g0_dot_keldysh(backward, a, w + 0.);
      }
    }

    chi = make_gf_from_fourier(chi_omega, time_mesh);
  };

  /// returns ......
  // alpha not supported yet
  mda::array<dcomplex, 1> operator()(gf_index_t const &a, gf_index_t const &b) const {
    if (a.spin == b.spin) {
      TRIQS_RUNTIME_ERROR << "a and b should have opposite spins.";
      //return 0.0 * model.omegas;
    }

    auto out = chi(b.contour.time - a.contour.time)(mda::range(), static_cast<int>(a.contour.k_idx),
                                                    static_cast<int>(b.contour.k_idx));
    return mda::exp(1.0i * omegas * (a.contour.time - time)) * out;
  };

  int get_nr_omegas() const { return omegas.size(); };
};

class integrand_spin_plus_spin_minus_freq {
  g0_keldysh_contour_t g0;
  gf_index_t g_idx_X; // Fixed Point
  chi_function_t const chi;

 public:
  /// Returns integrand for the specified times
  using result_t = array<dcomplex, 1>;
  [[nodiscard]] std::pair<result_t, int> operator()(std::vector<double> const &times,
                                                    bool const keep_u_hypercube = true) const;

  integrand_spin_plus_spin_minus_freq(g0_keldysh_contour_t g0_, double time, mda::array<double, 1> const &omegas)
     : g0(std::move(g0_)), g_idx_X{time, up, forward}, chi(g0.model.model_omega, omegas, time, backward){};
};

} // namespace keldy::impurity_oneband
