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

#include "keldy/warpers/plasma_uv.hpp"
#include "keldy/warpers/product_1d_simple.hpp"
#include "model.hpp"
#include "wick_direct.hpp"
#include "wick_kernel.hpp"

#include "../common.hpp"
#include "../integrator.hpp"
#include "../warpers/warpers.hpp"

#include "../interfaces/gsl_vegas_wrap.hpp"
#include "../interfaces/cuba_interface.hpp"

#include <triqs/utility/first_include.hpp>
#include <string>
#include <algorithm>

namespace keldy::impurity_oneband {

// ******************************************************************************************************************************************************
// Direct Evaluation ('Profumo')

inline warper_product_1d_simple_t simple_plasma_warper_factory(std::string const &label, integrand_g_direct const &f,
                                                               double time, int nr_sample_points_warper,
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
class compute_charge_Q_direct : public integrator<dcomplex, integrand_g_direct> {
 public:
  compute_charge_Q_direct(g0_model model, double time, int order, double cutoff_integrand,
                          std::string warper_function_name, int nr_sample_points_warper, double warper_scale = 1)
     : integrator{dcomplex{0},
                  integrand_g_direct{g0_keldysh_contour_t{model}, gf_index_t{time, up, forward},
                                     gf_index_t{time, up, backward}, cutoff_integrand},
                  {},
                  order,
                  "sobol",
                  0} {

    warper.warpers.emplace_back(warper_plasma_uv_t(time));
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
                  "sobol",
                  0} {

    warper.warpers.emplace_back(warper_plasma_uv_t(time));
    warper.warpers.emplace_back(warper_product_1d_t{fn_, time, nr_sample_points_warper});
  }
};

class CPP2PY_IGNORE adapt_integrand {
  // double time_max_;
  integrand_g_direct integrand_;
  warper_product_1d_simple_t pre_warper;

 public:
  adapt_integrand(double time_max, integrand_g_direct integrand, double warper_scale)
     : //time_max_(time_max),
       integrand_(std::move(integrand)),
       pre_warper{simple_plasma_warper_factory("inverse_square", integrand_, time_max, int(1e6), warper_scale)} {};

  double operator()(std::vector<double> const &li_vec) {
    std::vector<double> ui_vec = pre_warper.ui_from_li(li_vec);

    auto [eval, in_domain] = integrand_(ui_vec);
    eval *= pre_warper.jacobian(li_vec);

    int order = li_vec.size();
    return std::real(-std::pow(1_j, 1 + order) * eval);
  }
};

// Class to compute charge = G^{lesser}_{up,up}(t).
class compute_charge_Q_direct_gsl_vegas : public gsl_vegas_wrapper_t {
 public:
  compute_charge_Q_direct_gsl_vegas(model_param_t params, double time, int order, double cutoff_integrand,
                                    std::string gsl_rng_name, double warper_scale = 1)
     : gsl_vegas_wrapper_t{
        adapt_integrand{time,
                        integrand_g_direct{g0_keldysh_contour_t{g0_model{g0_model_omega{params}, false}},
                                           gf_index_t{time, up, forward}, gf_index_t{time, up, backward},
                                           cutoff_integrand},
                        warper_scale},
        order, 1.0, gsl_rng_name} {}
};

// Class to compute charge = G^{lesser}_{up,up}(t).
class compute_charge_Q_direct_cuba : public cuba_wrapper {
 public:
  compute_charge_Q_direct_cuba(model_param_t params, double time, int order, double cutoff_integrand,
                               cuba_common_param in, double warper_scale = 1)
     : cuba_wrapper{adapt_integrand{time,
                                    integrand_g_direct{g0_keldysh_contour_t{g0_model{g0_model_omega{params}, false}},
                                                       gf_index_t{time, up, forward}, gf_index_t{time, up, backward},
                                                       cutoff_integrand},
                                    warper_scale},
                    order, std::move(in)} {}
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
                  "sobol",
                  0} {
    warper.warpers.emplace_back(warper_plasma_uv_t(time));
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

inline std::function<double(double)> scalar_warper_function_factory_kernel(std::string const &label,
                                                                           integrand_g_kernel const &f, double time) {
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
                    int nr_sample_points_warper)
     : integrator{kernel_binner{0.0, time, 100},
                  integrand_g_kernel{g0_keldysh_contour_t{g0_model{g0_model_omega{params}, false}},
                                     gf_index_t{time, up, forward}},
                  {},
                  order,
                  "sobol",
                  0} {

    warper.warpers.emplace_back(warper_plasma_uv_t(time));
    warper.warpers.emplace_back(warper_product_1d_simple_t{
       scalar_warper_function_factory_kernel(warper_function_name, integrand, time), time, nr_sample_points_warper});
  }
};

} // namespace keldy::impurity_oneband
