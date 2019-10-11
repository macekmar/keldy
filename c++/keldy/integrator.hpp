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

template <typename W>
class CPP2PY_IGNORE integrator_t {
 public:
  W warper;

 private:
  std::function<void(std::vector<double> const &, double)> acc;
  std::function<std::vector<double>()> rng;
  mpi::communicator comm;

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
      acc(ui_vec, warper.jacobian(li_vec));
    }
  }

  integrator_t() = default;

  integrator_t(std::function<void(std::vector<double> const &, double)> acc_, W w, int dimension, const std::string& rng_name,
               int rng_seed, mpi::communicator comm_)
     : warper(std::move(w)), acc(std::move(acc_)), comm(std::move(comm_)) {
    if (rng_name == "sobol") {
      rng = sobol(dimension, rng_seed);
    } else {
      TRIQS_RUNTIME_ERROR << "No other rng available.";
    }
  }
};

} // namespace keldy
