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

#include "model.hpp"
#include "wick_direct.hpp"

#include "keldy/common.hpp"
#include "keldy/integrator.hpp"
#include "keldy/warpers/warpers.hpp"

#include <triqs/utility/first_include.hpp>

namespace keldy::impurity_oneband {

class compute_charge_Q {
 private:
  dcomplex result = 0.0;
  integrator_t<warper_plasma_simple_t> integrator;
  mpi::communicator comm;

 public:
  compute_charge_Q(int order, double time, model_param_t params, int nr_sample_points_ansatz) {

    g0_model g0_model(params);

    gf_index_t a(time, up, backward);
    gf_index_t b(time, up, forward);

    auto f = integrand_g_t1t2_direct{g0_keldysh_contour_t{g0_model}, a, b};

    auto f1 = [f](double t) { return std::abs(f(std::vector<double>{t})); };
    warper_plasma_simple_t warper{f1, params.time_max, nr_sample_points_ansatz};

    auto f2 = [&result = this->result, f](std::vector<double> const &ui_vec, double jac) { result += jac * f(ui_vec); };

    integrator = integrator_t{f2, warper, order, "sobol", comm};
  }

  void run(int nr_steps) { integrator.run(nr_steps); }

  // dcomplex reduce_result() const { return mpi::all_reduce(result, comm); }
  dcomplex reduce_result() const { return result; }

  // FIXME
  // warper_t get_warper() {return integrator.warper}

  // get_integrand (return std::function)
};

} // namespace keldy::impurity_oneband
