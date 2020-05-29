/*******************************************************************************
 *
 * keldy
 *
 * Copyright (c) 2019-2020 The Simons Foundation
 * Copyright (c) 2019-2020 CEA: Commissariat à l’énergie atomique
 *                              et aux énergies alternatives
 *   authors: Philipp Dumitrescu, Marjan Macek, Corentin Bertrand
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

#include "keldy/warpers/plasma_uv.hpp"
#include "keldy/warpers/product_1d_simple.hpp"
#include "model.hpp"
#include "wick_direct.hpp"
#include "wick_direct_time.hpp"
#include "wick_kernel.hpp"

#include "../common.hpp"
#include "../integrator.hpp"
#include "../warpers/warpers.hpp"

#include <triqs/utility/first_include.hpp>
#include <string>
#include <algorithm>
#include <utility>

namespace keldy::impurity_oneband {

// ******************************************************************************************************************************************************
// Direct Evaluation ('Profumo')

// Class to compute charge = G^{lesser}_{up,up}(t).
class compute_charge_Q_direct : public integrator<dcomplex, impurity_oneband::integrand_g_direct> {
 public:
  compute_charge_Q_direct(g0_model model, double time, int order, double cutoff_integrand)
     : integrator{dcomplex{0},
                  integrand_g_direct{g0_keldysh_contour_t{std::move(model)}, gf_index_t{time, up, forward},
                                     gf_index_t{time, up, backward}, cutoff_integrand},
                  {},
                  order,
                  "sobol_unshifted",
                  0} {}

  compute_charge_Q_direct(model_param_t params, double time, int order, double cutoff_integrand)
     : compute_charge_Q_direct{g0_model{g0_model_omega{params}, false}, time, order, cutoff_integrand} {};
};

// Class to compute the current = -2e/hbar gamma Re[G^<_{lead-dot}(t, t)]
class compute_current_J_direct : public integrator<dcomplex, integrand_g_direct> {
 public:
  compute_current_J_direct(g0_model model, double time, int order, double cutoff_integrand)
     : integrator{dcomplex{0},
                  integrand_g_direct{g0_keldysh_contour_t{std::move(model)}, gf_index_t{time, up, forward, 0, 1},
                                     gf_index_t{time, up, backward, 0, 0}, cutoff_integrand},
                  {},
                  order,
                  "sobol_unshifted",
                  0} {}

  compute_current_J_direct(model_param_t params, double time, int order, double cutoff_integrand)
     : compute_current_J_direct{g0_model{g0_model_omega{params}, true}, time, order, cutoff_integrand} {};
};
// ****************************************************************************
// Profumo + Time evolution methods

// Class to compute charge = G^{lesser}_{up,up}(t).
class compute_charge_Q_direct_time : public integrator<binner::binner_t<1>, integrand_g_direct_time> {
 public:
  compute_charge_Q_direct_time(g0_model model, double time, int order, int nr_time_slices, double cutoff_integrand)
     : integrator{binner::binner_t<1>({std::make_tuple(0., time, nr_time_slices)}),
                  integrand_g_direct_time{g0_keldysh_contour_t{std::move(model)}, gf_index_t{time, up, forward},
                                          gf_index_t{time, up, backward}, cutoff_integrand},
                  {},
                  order,
                  "sobol_unshifted",
                  0} {}

  compute_charge_Q_direct_time(model_param_t params, double time, int order, int nr_time_slices,
                               double cutoff_integrand)
     : compute_charge_Q_direct_time{g0_model{g0_model_omega{params}, false}, time, order, nr_time_slices,
                                    cutoff_integrand} {};
};

// ******************************************************************************************************************************************************
// Kernal Method

class compute_gf_kernel : public integrator<binner::binner_t<1, 1>, integrand_g_kernel> {
 public:
  compute_gf_kernel(g0_model model, double time, int order, int nr_bins = 100)
     : integrator{binner::binner_t<1, 1>({std::make_tuple(0.0, time, nr_bins)}, {2}),
                  integrand_g_kernel{g0_keldysh_contour_t{std::move(model)}, gf_index_t{time, up, forward}},
                  {},
                  order,
                  "sobol_unshifted",
                  0} {}

  compute_gf_kernel(model_param_t params, double time, int order, int nr_bins = 100)
     : compute_gf_kernel{g0_model{g0_model_omega{params}, false}, time, order, nr_bins} {}
};

} // namespace keldy::impurity_oneband
