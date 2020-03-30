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
#include "wick_kernel.hpp"

#include "../common.hpp"
#include "../integrator.hpp"
#include "../warpers/warpers.hpp"

#include <triqs/utility/first_include.hpp>
#include <string>
#include <algorithm>

namespace keldy::impurity_oneband {

// ******************************************************************************************************************************************************
// Direct Evaluation ('Profumo')

inline warpers::warper_product_1d_simple_t simple_plasma_warper_factory(std::string const &label,
                                                                        integrand_g_direct const &f, double time,
                                                                        int nr_sample_points_warper,
                                                                        double warper_scale) {
  if (label == "first_order") {
    return {[time, &f](double t) -> double { return std::abs(f(std::vector<double>{time - t}).first) + 1e-12; }, time,
            nr_sample_points_warper};
  }
  if (label == "inverse_square") {
    return {[warper_scale](double t) -> double {
              return warper_scale * warper_scale / ((warper_scale + t) * (warper_scale + t));
            },
            [warper_scale](double t) -> double { return warper_scale * t / (warper_scale + t); },
            [warper_scale](double l) -> double { return warper_scale * l / (warper_scale - l); }, time,
            nr_sample_points_warper};
  }
  if (label == "exponential") {
    return {[warper_scale](double t) -> double { return std::exp(-(t / warper_scale)); },
            [warper_scale](double t) -> double { return warper_scale * (1 - std::exp(-t / warper_scale)); },
            [warper_scale](double l) -> double { return -warper_scale * std::log(1 - l / warper_scale); }, time,
            nr_sample_points_warper};
  }
  if (label == "identity") {
    return {time};
  }
  TRIQS_RUNTIME_ERROR << "Warper function name is not defined.";
}

// Class to compute charge = G^{lesser}_{up,up}(t).
class compute_charge_Q_direct : public integrator<dcomplex, impurity_oneband::integrand_g_direct> {
 public:
  compute_charge_Q_direct(g0_model model, double time, int order, double cutoff_integrand,
                          std::string warper_function_name, int nr_sample_points_warper, double warper_scale = 1)
     : integrator{dcomplex{0},
                  integrand_g_direct{g0_keldysh_contour_t{model}, gf_index_t{time, up, forward},
                                     gf_index_t{time, up, backward}, cutoff_integrand},
                  {},
                  order,
                  "sobol_unshifted",
                  0} {

    warper.warpers.emplace_back(warpers::warper_plasma_uv_t(time));
    warper.warpers.emplace_back(
       simple_plasma_warper_factory(warper_function_name, integrand, time, nr_sample_points_warper, warper_scale));
  }

  compute_charge_Q_direct(model_param_t params, double time, int order, double cutoff_integrand,
                          std::string warper_function_name, int nr_sample_points_warper, double warper_scale = 1)
     : compute_charge_Q_direct{g0_model{g0_model_omega{params}, false},
                               time,
                               order,
                               cutoff_integrand,
                               warper_function_name,
                               nr_sample_points_warper,
                               warper_scale} {};
};

class compute_charge_Q_direct_plasma_1D : public integrator<dcomplex, integrand_g_direct> {
 public:
  compute_charge_Q_direct_plasma_1D(model_param_t params, double time, int order,
                                    std::vector<std::function<double(double)>> fn_, int nr_sample_points_warper)
     : integrator{dcomplex{0},
                  integrand_g_direct{g0_keldysh_contour_t{g0_model{g0_model_omega{params}, false}},
                                     gf_index_t{time, up, forward}, gf_index_t{time, up, backward}},
                  {},
                  order,
                  "sobol_unshifted",
                  0} {

    warper.warpers.emplace_back(warpers::warper_plasma_uv_t(time));
    warper.warpers.emplace_back(warpers::warper_product_1d_t{fn_, time, nr_sample_points_warper});
  }
};

// Class to compute the current = -2e/hbar gamma Re[G^<_{lead-dot}(t, t)]
class compute_current_J_direct : public integrator<dcomplex, integrand_g_direct> {
 public:
  compute_current_J_direct(g0_model model, double time, int order, double cutoff_integrand,
                           std::string warper_function_name, int nr_sample_points_warper, double warper_scale = 1)
     : integrator{dcomplex{0},
                  integrand_g_direct{g0_keldysh_contour_t{model}, gf_index_t{time, up, forward, 0, 1},
                                     gf_index_t{time, up, backward, 0, 0}, cutoff_integrand},
                  {},
                  order,
                  "sobol_unshifted",
                  0} {
    warper.warpers.emplace_back(warpers::warper_plasma_uv_t(time));
    warper.warpers.emplace_back(
       simple_plasma_warper_factory(warper_function_name, integrand, time, nr_sample_points_warper, warper_scale));
  }

  compute_current_J_direct(model_param_t params, double time, int order, double cutoff_integrand,
                           std::string warper_function_name, int nr_sample_points_warper, double warper_scale = 1)
     : compute_current_J_direct{g0_model{g0_model_omega{params}, true},
                                time,
                                order,
                                cutoff_integrand,
                                warper_function_name,
                                nr_sample_points_warper,
                                warper_scale} {};
};
// ******************************************************************************************************************************************************
// Kernal Method

inline warpers::warper_product_1d_simple_t simple_plasma_warper_factory_kernel(std::string const &label,
                                                                               integrand_g_kernel const &f, double time,
                                                                               int nr_sample_points_warper,
                                                                               double warper_scale) {
  if (label == "inverse_square") {
    return {[warper_scale](double t) -> double {
              return warper_scale * warper_scale / ((warper_scale + t) * (warper_scale + t));
            },
            [warper_scale](double t) -> double { return warper_scale * t / (warper_scale + t); },
            [warper_scale](double l) -> double { return warper_scale * l / (warper_scale - l); }, time,
            nr_sample_points_warper};
  }
  if (label == "exponential") {
    return {[warper_scale](double t) -> double { return std::exp(-(t / warper_scale)); },
            [warper_scale](double t) -> double { return warper_scale * (1 - std::exp(-t / warper_scale)); },
            [warper_scale](double l) -> double { return -warper_scale * std::log(1 - l / warper_scale); }, time,
            nr_sample_points_warper};
  }
  if (label == "identity") {
    return {time};
  }
  TRIQS_RUNTIME_ERROR << "Warper function name is not defined.";
}

inline std::function<double(double)> scalar_warper_function_factory_kernel(std::string const &label,
                                                                           integrand_g_kernel const &f,
                                                                           double time) CPP2PY_IGNORE {
  if (label == "first_order") {
    return [time, &f](double t) -> double {
      return f(std::vector<double>{time - t}).first.sum_weights() + 1e-12;
    }; // refine
  }
  if (label == "identity") {
    return []([[maybe_unused]] double t) -> double { return 1.0; };
  }
  TRIQS_RUNTIME_ERROR << "Warper function name is not defined.";
}

class compute_gf_kernel : public integrator<kernel_binner, integrand_g_kernel> {
 public:
  compute_gf_kernel(model_param_t params, double time, int order, std::string warper_function_name,
                    int nr_sample_points_warper, int nr_bins = 100)
     : integrator{kernel_binner{0.0, time, nr_bins},
                  integrand_g_kernel{g0_keldysh_contour_t{g0_model{g0_model_omega{params}, false}},
                                     gf_index_t{time, up, forward}},
                  {},
                  order,
                  "sobol_unshifted",
                  0} {

    warper.warpers.emplace_back(warpers::warper_plasma_uv_t(time));
    warper.warpers.emplace_back(warpers::warper_product_1d_simple_t{
       scalar_warper_function_factory_kernel(warper_function_name, integrand, time), time, nr_sample_points_warper});
  }

  compute_gf_kernel(g0_model model, double time, int order, std::string warper_function_name,
                    int nr_sample_points_warper, double warper_scale, int nr_bins = 100)
     : integrator{kernel_binner{0.0, time, nr_bins},
                  integrand_g_kernel{g0_keldysh_contour_t{model}, gf_index_t{time, up, forward}},
                  {},
                  order,
                  "sobol_unshifted",
                  0} {

    warper.warpers.emplace_back(warpers::warper_plasma_uv_t(time));
    warper.warpers.emplace_back(simple_plasma_warper_factory_kernel(warper_function_name, integrand, time,
                                                                    nr_sample_points_warper, warper_scale));
  }
};

} // namespace keldy::impurity_oneband
