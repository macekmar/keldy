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

#include "gsl_vegas_wrap.hpp"

namespace keldy::extern_c {

extern "C" {

double gsl_vegas_wrapper_t_f_wrap(double x[], size_t dim, void *p) {
  auto this_ptr = (gsl_vegas_wrapper_t *)p;
  return this_ptr->f_(std::vector(x, x + dim));
}
}

} // namespace keldy::extern_c
