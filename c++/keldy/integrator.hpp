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

namespace keldy {

template <typename T>
constexpr bool is_binned_variable = false;
// if (is_binned_variable<R>) {}

template <typename R, typename I, typename W>
class integrator {
 protected:
  mpi::communicator comm{};

  // Result of the integration and number of samples
  R result = 0; // TODO: do = 0 this for kernel
  uint64_t n_points = 0;

  // Function that when called with a vector of times returns the integrand
  I integrand;

  // Warper which defines the sample points transofrm [li -> ui] & Jacobian
  W warper;

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

      result += warper.jacobian(li_vec) * integrand(ui_vec);
      n_points++;
    }
  }

  integrator(I integrand_, W warper_, int dimension, const std::string &rng_name, int rng_seed)
     : integrand(std::move(integrand_)), warper(std::move(warper_)) {
    if (rng_name == "sobol") {
      rng = sobol(dimension, rng_seed);
    } else {
      TRIQS_RUNTIME_ERROR << "No other rng available.";
    }
  }

  R reduce_result() const {
    R result_all = mpi::all_reduce(result, comm);
    return result_all / reduce_nr_points_run();
  }

  uint64_t reduce_nr_points_run() const { return mpi::all_reduce(n_points, comm); }

  // For Python Visualzation Purposes:
  auto get_integrand() const { return integrand; }
  auto get_warper() const { return warper; }
};

} // namespace keldy
