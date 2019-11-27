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
  int max_number_evaluations = int(1e10);
};

struct CPP2PY_IGNORE cuba_vegas_param {
  int n_start = 1000;
  int n_increase = 1000;
  int n_batch = 500;
  int gridno = 0;
};

CPP2PY_ARG_AS_DICT inline void fake1(cuba_output const &temp){};
CPP2PY_ARG_AS_DICT inline void fake2(cuba_common_param const &temp){};
CPP2PY_ARG_AS_DICT inline void fake3(cuba_vegas_param const &temp){};

// Do only for scalar functions

class cuba_vegas_wrapper {
 public:
  std::function<double(std::vector<double>)> f; // must make public for c interface hack

 private:
  cuba_common_param in;
  cuba_vegas_param in_v;
  cuba_output out{};

 public:
  cuba_vegas_wrapper(std::function<double(std::vector<double>)> f_, int dim, cuba_common_param in_, cuba_vegas_param in_v_)
     : f(std::move(f_)), in(std::move(in_)), in_v(std::move(in_v_)) {

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

  void run() {
    Vegas(in.n_dim, in.n_components, extern_c::cuba_f_wrap, this, in.n_points_vectorization, in.error_eps_rel,
          in.error_eps_abs, in.flags, in.seed, in.min_number_evaluations, in.max_number_evaluations, in_v.n_start,
          in_v.n_increase, in_v.n_batch, in_v.gridno, nullptr, nullptr, &out.n_eval, &out.error_flag, &out.result,
          &out.error, &out.chi_sq_prob);
  }

  cuba_output get_output() { return out; }
};

} // namespace keldy
