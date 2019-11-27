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

#include <complex>

namespace keldy {

// Parameters / definitions used in both real and imag time calculations
using dcomplex = std::complex<double>;
using time_real_t = double;
using time_imag_t = double;
using orbital_t = int;

enum spin_t { up = 0, down = 1};
enum keldysh_idx_t { forward = 0, backward = 1};


} // namespace keldy
