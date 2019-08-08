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
#include <mpi/mpi.hpp>
#include <omp.h>

namespace keldy {

template <typename W, typename M>
class CPP2PY_IGNORE integrator_t {
  W warper;

  std::function<void(M &, std::vector<double> const &, double)> acc;
  std::function<std::vector<double>()> rng;
  mpi::communicator comm;

 public:
  uint64_t run(M &measure, int nr_steps) {
    uint64_t n_pts = 0;
// No "for" in pragma: all threads execute full loop.
// Reduction adds to initial value of measure at the end : TODO
#pragma omp parallel reduction(+ : measure, n_pts)
    {
      auto local_rng = rng;
      int thread_num = omp_get_thread_num(), number_of_threads = omp_get_num_threads(); // call only once
      for (int i = 0; i < nr_steps; i++) {
        auto li_vec = local_rng();
        // if (i % comm.size() != comm.rank()) continue;
        if (i % number_of_threads != thread_num) continue;
        std::vector<double> ui_vec = warper.ui_from_li(li_vec);
        acc(measure, ui_vec, warper.jacobian(li_vec));
        ++n_pts;
      }
#pragma omp master
      { rng = local_rng; }
    }
    return n_pts;
  }

  uint64_t run_mpi(M &measure, int nr_steps) {
    int mpi_rank = comm.rank(), mpi_size = comm.size(); // call only once
    for (int i = 0; i < nr_steps; i++) {
      auto li_vec = rng();
      if (i % mpi_size != mpi_rank) continue;
      std::vector<double> ui_vec = warper.ui_from_li(li_vec);
      acc(measure, ui_vec, warper.jacobian(li_vec));
    }
    return nr_steps;
  }

  integrator_t() = default;

  integrator_t(std::function<void(M &, std::vector<double> const &, double)> acc_, W w, int dimension,
               std::string rng_name, mpi::communicator comm_)
     : warper(std::move(w)), acc(std::move(acc_)), comm(std::move(comm_)) {
    if (rng_name == "sobol") {
      rng = sobol(dimension);
    } else {
      TRIQS_RUNTIME_ERROR << "No other rng available.";
    }
  }
};

} // namespace keldy
