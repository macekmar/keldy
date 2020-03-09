//******************************************************************************
//
// keldy
//
// Copyright (C) 2019, The Simons Foundation
// authors: Philipp Dumitrescu, Olivier Parcollet
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

#include "common.hpp"
#include "qrng.hpp"
#include "warpers/warpers.hpp"
#include <cstddef>
#include <mpi/mpi.hpp>
#include <omp.h>
#include <triqs/utility/macros.hpp>
#include <triqs/arrays.hpp>

namespace keldy {

template <typename T>
constexpr bool is_binned_variable = false;
// if (is_binned_variable<R>) {}

template <typename R, typename I>
class integrator {
 protected:
  mpi::communicator comm{};

  int dimension = 1;

  // Result of the integration and number of samples
  R result{};
  uint64_t n_points = 0;
  uint64_t n_integrand_evals_in_domain = 0;

  // Function that when called with a vector of times returns the integrand
  I integrand;

  // Warpers which define the transform of sample points [li -> ui] & Jacobian
  warper_train_t warper{};

  // Random number generator
  std::function<std::vector<double>()> rng;

  // Checks that result of integrand can be added to result
  static_assert(std::is_same_v<decltype(std::declval<R &>() += std::declval<typename I::result_t>()), R &>);

  // TODO: static_assert that R can be mpi::all_reduce

 public:
  void run(int nr_steps) {
    int mpi_rank = comm.rank();
    int mpi_size = comm.size();
    for (int i = 0; i < nr_steps; i++) {
      auto li_vec = rng();
      if (i % mpi_size != mpi_rank) {
        continue;
      }
      std::vector<double> ui_vec = warper.ui_from_li(li_vec);

      auto [eval, in_domain] = integrand(ui_vec);
      eval *= warper.jacobian(li_vec);

      n_integrand_evals_in_domain += in_domain;
      result += eval;
      n_points++;
    }
  }

  triqs::arrays::array<double, 2> gather_data(int nr_points, std::function<std::vector<double>()> rng) {
    dcomplex im_unit (0.0,1.0);
    triqs::arrays::array<double, 2> values(nr_points, 5);
    values() = 0;

    int mpi_rank = comm.rank();
    int mpi_size = comm.size();
    
    for (int i = 0; i < nr_points; i++) {
      auto li_vec = rng();
      if (i % mpi_size != mpi_rank) {
        continue;
      }
      std::vector<double> ui_vec = warper.ui_from_li(li_vec);

      auto [eval, in_domain] = integrand(ui_vec);
      auto [min, max] = std::minmax_element(ui_vec.begin(), ui_vec.end());
      values(i,0) = *min;
      values(i,1) = *max;
      values(i,2) = ((-im_unit)*std::pow(im_unit, ui_vec.size())*eval).real();
      values(i,3) = in_domain;
      values(i,4) = warper.jacobian(li_vec);
    }

    values = mpi::mpi_all_reduce(values, comm);


    return values;
  }

  void reset_rng(std::string rng_name, int rng_state_seed, bool do_shift = false, bool do_scramble = false,
                 int rng_seed_shift = 0) {
    if (rng_name == "sobol_unshifted") {
      rng = sobol(dimension, rng_state_seed);
    } else if (rng_name == "sobol") {
      rng = sobol(dimension, rng_state_seed, do_shift, do_scramble, rng_seed_shift);
    } else if (rng_name == "pseudo") {
      rng = pseudo(dimension, rng_state_seed);
    } else {
      TRIQS_RUNTIME_ERROR << "No other rng available.";
    }
  }

  integrator(R result_, I integrand_, warper_train_t warper_, int dimension_, std::string rng_name, int rng_state_seed)
     : dimension(dimension_), result(std::move(result_)), integrand(std::move(integrand_)), warper(std::move(warper_)) {
    if (dimension <= 0) {
      TRIQS_RUNTIME_ERROR << "Dimension of the integration space must be > 0.";
    }
    reset_rng(rng_name, rng_state_seed);
  }

  R reduce_result() const {
    R result_all = mpi::all_reduce(result, comm);
    result_all *= 1.0 / reduce_nr_points_run();
    return result_all;
  }

  uint64_t reduce_nr_points_run() const { return mpi::all_reduce(n_points, comm); }

  uint64_t reduce_nr_points_in_domain() const { return mpi::all_reduce(n_integrand_evals_in_domain, comm); }

  // For Python Visualzation Purposes:
  auto get_integrand() const { return integrand; }
  auto get_warper() const { return warper; }
};

} // namespace keldy
