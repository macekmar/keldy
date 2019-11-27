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

inline std::function<double(double)>
scalar_warper_function_factory(std::string const &label, integrand_g_direct const &f, double time) CPP2PY_IGNORE {
  if (label == "first_order") {
    return [time, &f](double t) -> double { return std::abs(f(std::vector<double>{time - t})) + 1e-12; };
  }
  if (label == "lorentzian") {
    return [](double t) -> double { return 1.0 / ((1.0 + t) * (1.0 + t)); };
  }
  if (label == "identity") {
    return [](double t) -> double { return 1.0; };
  }
  TRIQS_RUNTIME_ERROR << "Warper function name is not defined.";
}

// Class to compute charge = G^{lesser}_{up,up}(t).
class compute_charge_Q_direct : public integrator<dcomplex, integrand_g_direct, warper_plasma_simple_t> {
 public:
  compute_charge_Q_direct(model_param_t params, double time, int order, std::string warper_function_name,
                          int nr_sample_points_warper)
     : integrator{dcomplex{0},
                  integrand_g_direct{g0_keldysh_contour_t{g0_model{params}}, gf_index_t{time, up, forward},
                                     gf_index_t{time, up, backward}},
                  warper_plasma_simple_t{time},
                  order,
                  "sobol",
                  0} {
    warper = {scalar_warper_function_factory(warper_function_name, integrand, time), time, nr_sample_points_warper};
  }
};

class CPP2PY_IGNORE adapt_integrand {
  double time_max_;
  integrand_g_direct integrand_;

 public:
  adapt_integrand(double time_max, integrand_g_direct integrand)
     : time_max_(time_max), integrand_(std::move(integrand)){};

  double operator()(std::vector<double> const &times_norm) {
    if (!std::is_sorted(std::begin(times_norm), std::end(times_norm))) {
      return 0.0;
    }
    auto times = times_norm;
    for (auto &t : times) {
      t *= time_max_;
    }
    int order = times.size();
    return std::real(-std::pow(1_j, 1 + order) * integrand_(times) * std::pow(time_max_, order));
  }
};

// Class to compute charge = G^{lesser}_{up,up}(t).
class compute_charge_Q_direct_gsl_vegas : public gsl_vegas_wrapper_t {
 public:
  compute_charge_Q_direct_gsl_vegas(model_param_t params, double time, int order, std::string gsl_rng_name)
     : gsl_vegas_wrapper_t{
        adapt_integrand{time,
                        integrand_g_direct{g0_keldysh_contour_t{g0_model{params}}, gf_index_t{time, up, forward},
                                           gf_index_t{time, up, backward}}},
        order, 1.0, gsl_rng_name} {}
};

// Class to compute charge = G^{lesser}_{up,up}(t).
class compute_charge_Q_direct_cuba_vegas : public cuba_vegas_wrapper {
 public:
  compute_charge_Q_direct_cuba_vegas(model_param_t params, double time, int order, cuba_common_param in,
                                     cuba_vegas_param in_v)
     : cuba_vegas_wrapper{
        adapt_integrand{time,
                        integrand_g_direct{g0_keldysh_contour_t{g0_model{params}}, gf_index_t{time, up, forward},
                                           gf_index_t{time, up, backward}}},
        order, std::move(in), std::move(in_v)} {}
};

// ******************************************************************************************************************************************************
// Kernal Method

inline std::function<double(double)> scalar_warper_function_factory_kernel(std::string const &label,
                                                                           integrand_g_kernel const &f,
                                                                           double time) CPP2PY_IGNORE {
  if (label == "first_order") {
    return [time, &f](double t) -> double { return f(std::vector<double>{time - t}).sum_weights() + 1e-12; }; // refine
  }
  if (label == "identity") {
    return [](double t) -> double { return 1.0; };
  }
  TRIQS_RUNTIME_ERROR << "Warper function name is not defined.";
}

class compute_gf_kernel : public integrator<kernel_binner, integrand_g_kernel, warper_plasma_simple_t> {
 public:
  compute_gf_kernel(model_param_t params, double time, int order, std::string warper_function_name,
                    int nr_sample_points_warper)
     : integrator{kernel_binner{0.0, time, 100},
                  integrand_g_kernel{g0_keldysh_contour_t{g0_model{params}}, gf_index_t{time, up, forward}},
                  warper_plasma_simple_t{time},
                  order,
                  "sobol",
                  0} {
    warper = {scalar_warper_function_factory_kernel(warper_function_name, integrand, time), time,
              nr_sample_points_warper};
  }
};

} // namespace keldy::impurity_oneband
