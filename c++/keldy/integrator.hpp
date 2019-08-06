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

#include "common.hpp"
#include "qrng.hpp"
#include "warpers/warpers.hpp"
#include <mpi/mpi.hpp>
#include <omp.h>

namespace keldy {

template <typename W, typename M>
class CPP2PY_IGNORE integrator_t {
  W warper; //make erased warper_t ?

  // f = accumulate(times_vector, jaobian)
  std::function<void(M&, std::vector<double> const &, double)> acc;
  std::function<std::vector<double>()> rng;
  mpi::communicator comm;

 public:

  void run(M& measure, int nr_steps) {
    // No "for" in pragma: all threads execute full loop. 
    // Reduction adds to initial value of measure at the end : TODO
    #pragma omp parallel reduction(+: measure)
    {
      auto local_rng = rng;
      for (int i = 0; i < nr_steps; i++) {
        auto li_vec = local_rng();
        // if (i % comm.size() != comm.rank()) continue;
        if (i % omp_get_num_threads() != omp_get_thread_num()) continue;
        std::vector<double> ui_vec = warper.ui_from_li(li_vec);
        acc(measure, ui_vec, warper.jacobian(li_vec));
      }
      #pragma omp master 
      {
        rng = local_rng;
      }
    }
  }

  integrator_t(){};

  integrator_t(std::function<void(M&, std::vector<double> const &, double)> acc_,
               W w, int dimension, std::string rng_name, mpi::communicator comm_)
     : warper(std::move(w)), acc(std::move(acc_)), comm(std::move(comm_)) {
    if (rng_name == "sobol") {
      rng = sobol(dimension);
    } else {
      TRIQS_RUNTIME_ERROR << "No other rng available.";
    }
  }
};

} // namespace keldy
