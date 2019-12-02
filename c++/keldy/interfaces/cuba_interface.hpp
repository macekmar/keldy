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

#include <cuba.h>
#include <triqs/utility/exceptions.hpp>
#include <triqs/utility/first_include.hpp>
#include <vector>
#include <bitset>

namespace keldy {

namespace extern_c {
extern "C" {
int cuba_f_wrap(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata);
}
} // namespace extern_c

struct CPP2PY_IGNORE cuba_output {
  int n_regions = 0;
  int n_eval = 0;
  int error_flag = 0;

  double result = 0.0;
  double error = 0.0;
  double chi_sq_prob = 0.0;
};

struct CPP2PY_IGNORE cuba_common_param {
  int n_dim; // will be overwritten
  int n_components = 1;
  int n_points_vectorization = 1;
  double error_eps_rel = 1e-12;
  double error_eps_abs = 1e-12;

  int flags = 0;                              // Gets overwritten
  int verbosity = 0;                          // 0-3
  bool use_last_sampleset_only = false;       // 0 = off , 1 = on
  bool sample_function_smoothing_off = false; // 0 =  do smoothing , 1 = no smoothing
  bool store_state_after_run = false;         // 0 off / 1 store state

  std::string rng_type = "sobol"; // "sobol", "mt", "ranlux"
  int seed = 1;
  int randlux_level = 0;

  int min_number_evaluations = 1000;
  int max_number_evaluations = int(1e9);
};

struct CPP2PY_IGNORE cuba_vegas_param {
  int n_evals_per_iteration_start = 1000;
  int n_evals_per_iteration_increase = 500;
  int n_samples_per_batch = 1000;
  int internal_store_grid_nr = 0;
};

struct CPP2PY_IGNORE cuba_suave_param {
  int n_new_evals_each_subdivision = 5000;
  int n_min_samples_region_threashold = 500;
  double flatness_parameter_p = 50.0;
};

// struct CPP2PY_IGNORE cuba_divonne_param {
//   int key_1_sampling_during_partition = -1;
//   int key_2_sampling_during_final_integration = -1;
//   int key_3_strategy_refinement = 1;
//   int max_pass_partitioning = 5;
//   double border_width = 0.0;
//   double max_chi_sq_per_subregion = 5.0;
//   double min_deviation = 1.0;
//   // int n_points_given_array;
//   // iter
//   // weight
// };

CPP2PY_ARG_AS_DICT inline void fake_output(cuba_output const &temp){};
CPP2PY_ARG_AS_DICT inline void fake_common(cuba_common_param const &temp){};
CPP2PY_ARG_AS_DICT inline void fake_vegas(cuba_vegas_param const &temp){};
CPP2PY_ARG_AS_DICT inline void fake_suave(cuba_suave_param const &temp){};
// CPP2PY_ARG_AS_DICT inline void fake_divonne(cuba_divonne_param const &temp){};

// Do only for scalar functions

class cuba_wrapper {
 public:
  std::function<double(std::vector<double>)> f; // must make public for c interface hack

 private:
  cuba_common_param in;
  cuba_output out{};

 public:
  cuba_wrapper(std::function<double(std::vector<double>)> f_, int dim, cuba_common_param in_)
     : f(std::move(f_)), in(std::move(in_)) {

    in.n_dim = dim;

    if (in.rng_type == "sobol") {
      in.seed = 0;
      in.randlux_level = 0;
    } else if (in.rng_type == "mt") {
      if (in.seed == 0) {
        TRIQS_RUNTIME_ERROR << "seed can't be 0 for mt";
      }
      in.randlux_level = 0;
    } else if (in.rng_type == "ranlux") {
      if (in.seed == 0) {
        TRIQS_RUNTIME_ERROR << "seed can't be 0 for ranlux";
      }
      if (in.randlux_level < 0) {
        TRIQS_RUNTIME_ERROR << "randlux_level should be > 0";
      }
      if (in.randlux_level == 0) {
        in.randlux_level = 24;
      }
    } else {
      TRIQS_RUNTIME_ERROR << "rng_type unknown";
    }

    if (in.verbosity < 0 || in.verbosity > 3) {
      TRIQS_RUNTIME_ERROR << "verbosity Out of Bounds";
    }

    in.flags = (in.verbosity << 0) + (int(in.use_last_sampleset_only) << 2)
       + (int(in.sample_function_smoothing_off) << 3) + (int(in.store_state_after_run) << 4) + (in.randlux_level << 8);
  }

  void run_vegas(cuba_vegas_param in_v) {
    Vegas(in.n_dim, in.n_components, extern_c::cuba_f_wrap, this, in.n_points_vectorization, in.error_eps_rel,
          in.error_eps_abs, in.flags, in.seed, in.min_number_evaluations, in.max_number_evaluations,
          in_v.n_evals_per_iteration_start, in_v.n_evals_per_iteration_increase, in_v.n_samples_per_batch,
          in_v.internal_store_grid_nr, nullptr, nullptr, &out.n_eval, &out.error_flag, &out.result, &out.error,
          &out.chi_sq_prob);
  }

  void run_suave(cuba_suave_param in_s) {
    Suave(in.n_dim, in.n_components, extern_c::cuba_f_wrap, this, in.n_points_vectorization, in.error_eps_rel,
          in.error_eps_abs, in.flags, in.seed, in.min_number_evaluations, in.max_number_evaluations,
          in_s.n_new_evals_each_subdivision, in_s.n_min_samples_region_threashold, in_s.flatness_parameter_p, nullptr,
          nullptr, &out.n_regions, &out.n_eval, &out.error_flag, &out.result, &out.error, &out.chi_sq_prob);
  }

  // TODO: Wrapper
  // void run_divonne() {
  //   Divonne(in.n_dim, in.n_components, extern_c::cuba_f_wrap, this, in.n_points_vectorization, in.error_eps_rel,
  //           in.error_eps_abs, in.flags, in.seed, in.min_number_evaluations, in.max_number_evaluations, const int key1,
  //           const int key2, const int key3, const int maxpass, const double border, const double maxchisq,
  //           const double mindeviation, const int ngiven, const int ldxgiven, double xgiven[], const int nextra,
  //           peakfinder_t peakfinder, const char *statefile, void *spin, &out.n_regions, &out.n_eval, &out.error_flag,
  //           &out.result, &out.error, &out.chi_sq_prob);
  // }


  // key_integration_order: 
  //      k = 7,9,11,13 else defaults. Degree of cubature rule. 
  //      k = 13 only for dim = 2, k = 11 only for dim = 3
  //      defaults to max key available for given dim
  void run_cuhre(int key_integration_order) {
    Cuhre(in.n_dim, in.n_components, extern_c::cuba_f_wrap, this, in.n_points_vectorization, in.error_eps_rel,
          in.error_eps_abs, in.flags, in.min_number_evaluations, in.max_number_evaluations, key_integration_order,
          nullptr, nullptr, &out.n_regions, &out.n_eval, &out.error_flag, &out.result, &out.error, &out.chi_sq_prob);
  }

  cuba_output get_output() { return out; }
};

} // namespace keldy
