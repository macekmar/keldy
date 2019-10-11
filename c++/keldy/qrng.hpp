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
#include "qrng_details/digitalseq_b2g.hpp"
#include <triqs/utility/exceptions.hpp>

namespace keldy {

// // rng type interface -- copied from std::mersenne_twister_engine
// class rng_t {
//  public:
//   std::vector<double> operator()(); // make buffered function... (need skip when threading with MPI) in triqs
//   void seed();
//   void discard();
//   // ctor
// };

class new_class {};

class sobol {
 public:
  std::vector<double> operator()() {
    if (generator.past_end()) {
      TRIQS_RUNTIME_ERROR << "More sample points requested then generator can hold.";
    }
    // figure out MPI threading (i.e. discard(n))
    std::vector<double> result = generator.next();
    return result;
  }

  // k = 0 is like reset
  void seed(int k) { generator.set_state(std::uint32_t(k)); }

  void discard(int nr_discard) {
    if (nr_discard != 0) {
      generator.set_state(generator.k + nr_discard);
    }
  }

  sobol(int dim, int rng_seed, int log_max_points_ = 31)
     : //dimension(dim),
       log_max_points(log_max_points_),
       generator(qmc::JK2008_sobolseq<double, std::uint32_t>(dim, log_max_points_)) {
    if (log_max_points > 31) {
      TRIQS_RUNTIME_ERROR << "To many points requested from Sobol generator.";
    }
    generator.set_state(std::uint32_t(rng_seed));
  }

 private:
  //int dimension;
  int log_max_points;
  qmc::digitalseq_b2g<double, std::uint32_t> generator;
};

// Harmonic Generator

} // namespace keldy
