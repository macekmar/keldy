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

#include <string>
#include <gsl/gsl_rng.h>
#include <triqs/utility/first_include.hpp>

namespace keldy {

class CPP2PY_IGNORE gsl_rng_wrapper_t {
  gsl_rng *rng;

 public:
  gsl_rng_wrapper_t(std::string const &gsl_rng_name) {
    if (gsl_rng_name == "mt") {
      rng = gsl_rng_alloc(gsl_rng_mt19937);
    }
    if (gsl_rng_name == "ranlxs0") {
      rng = gsl_rng_alloc(gsl_rng_ranlxs0);
    }
    if (gsl_rng_name == "cmrg") {
      rng = gsl_rng_alloc(gsl_rng_cmrg);
    }
    if (gsl_rng_name == "taus2") {
      rng = gsl_rng_alloc(gsl_rng_taus2);
    }

    // Read Default Environment Set-Up
    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_default);
  }

  void seed(unsigned long int s) { gsl_rng_set(rng, s); }

  auto get_ptr() { return rng; }

  ~gsl_rng_wrapper_t() { gsl_rng_free(rng); }
};

} // namespace keldy
