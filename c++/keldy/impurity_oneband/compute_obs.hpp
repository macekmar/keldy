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

#include "keldy/common.hpp"
#include "keldy/integrator.hpp"
#include "keldy/warpers/warpers.hpp"

#include <triqs/utility/first_include.hpp>

namespace keldy::impurity_oneband {

class CPP2PY_IGNORE output_scalar_t {
 public:
  dcomplex result = 0.0;
  int nr_points = 0;



  // for reduction
  output_scalar_t &operator+=(const output_scalar_t &rhs) {
    this->result += rhs.result;
    this->nr_points += rhs.nr_points;
    return *this;
  }

  // output_scalar_t() = default;
};

#pragma omp declare reduction(+ : output_scalar_t : omp_out += omp_in) \
  initializer(omp_priv = output_scalar_t{})


class compute_charge_Q {

 private:
  output_scalar_t output;
  integrator_t<warper_plasma_simple_t, output_scalar_t> integrator;
  mpi::communicator comm;

 public:
  integrand_g_t1t2_direct integrand; // keep public copy for viz purposes

  compute_charge_Q(int order, double time, model_param_t params, int nr_sample_points_ansatz)
     : integrand{g0_keldysh_contour_t{g0_model{params}}, gf_index_t{time, up, backward},
                 gf_index_t{time, up, forward}} {

    auto f1 = [time, f = this->integrand](double t) { return std::abs(f(std::vector<double>{time - t})) + 1e-12; };

    warper_plasma_simple_t warper{f1, time, nr_sample_points_ansatz};

    auto f2 = [f = this->integrand](output_scalar_t &out, std::vector<double> const &ui_vec, double jac) {
      out.result += jac * f(ui_vec);
      out.nr_points++;
    };

    integrator = integrator_t<warper_plasma_simple_t, output_scalar_t>{f2, warper, order, "sobol", comm};
  }

  void run(int nr_steps) { 
    integrator.run(output, nr_steps); 
  }

  dcomplex reduce_result() const {
    dcomplex result_all = mpi::all_reduce(output.result, comm);
    return result_all / get_nr_points_run();
  }

  int get_nr_points_run() const {
    int nr_points_total = mpi::all_reduce(output.nr_points, comm);
    return nr_points_total;
  }

  // FIXME
  // warper_t get_warper() {return integrator.warper}
};

} // namespace keldy::impurity_oneband
