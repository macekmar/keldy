/*******************************************************************************
 *
 * keldy
 *
 * Copyright (c) 2020 The Simons Foundation
 * Copyright (c) 2020 CEA: Commissariat à l’énergie atomique
 *                         et aux énergies alternatives
 *   authors: Philipp Dumitrescu, Marjan Macek
 *
 * keldy is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * keldy is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * keldy. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include "warpers_common.hpp"
#include "product_1d_simple.hpp"
#include "product_1d.hpp"
#include "../interfaces/gsl_fit_linear_wrap.hpp"
#include "../interfaces/gsl_interp_wrap.hpp"
#include "../interfaces/gsl_minimize.hpp"
#include "../qrng.hpp"
#include "../binner.hpp"

#include <triqs/arrays.hpp>
#include <triqs/arrays/mpi.hpp>
//#include <mpi/mpi.hpp>

#include <algorithm>
#include <functional>

namespace keldy::warpers {

namespace nda = triqs::arrays;

namespace proj_details {

auto gaussian = [](double x, double x0, double sigma) { return std::exp(-std::pow((x - x0) / sigma, 2) / 2); };

// sigma is in unit of time step
inline array<double, 1> CPP2PY_IGNORE kernel_smoothing(array<double, 1> const &y_in, double const sigma) {
  int const N = y_in.size();
  array<double, 1> log_y_out(N);
  array<double, 1> x(N);
  std::vector<int> dead_points;
  for (int i = 0; i < N; ++i) {
    x(i) = i;
    if (y_in(i) <= 0) {
      dead_points.emplace_back(i);
    }
  }

  int const H = std::max(std::min(int(3 * sigma), N), 1);
  int const K = 2 * H + 1;
  array<double, 1> kernel(K);
  for (int j = 0; j < K; ++j) {
    kernel(j) = gaussian(j, H, sigma);
  }

  details::local_linear_reg(x, log(y_in), log_y_out, kernel, dead_points);

  //for (int k = 0; k < N; ++k) {
  //  if (y_out(k) < 0) y_out(k) = 0; // TODO: or should we raise an error?
  //}

  return exp(log_y_out);
}

void CPP2PY_IGNORE gather_data(std::function<dcomplex(std::vector<double>)> const &warped_integrand,
                               std::vector<binner::binner_t<1, 0, double>> &xi, int const nr_samples, int order) {
  double value;
  mpi::communicator comm{};
  int mpi_rank = comm.rank();
  int mpi_size = comm.size();

  sobol g = sobol(order, 0); // TODO: different generators
  for (int i = 0; i < nr_samples; i++) {
    auto l = g();

    if (i % mpi_size != mpi_rank) {
      continue;
    }

    value = std::abs(warped_integrand(l));

    for (int axis = 0; axis < order; axis++) {
      xi[axis](l[axis]) << value;
    }
  }
  for (int axis = 0; axis < order; axis++) {
    xi[axis] = mpi_reduce(xi[axis], comm, 0, true);
  }

  for (int axis = 0; axis < order; axis++) {
    auto data = xi[axis].data();             // view
    auto count = xi[axis].nr_values_added(); // view
    for (int bin = 0; bin < xi[axis].get_nr_bins(); bin++) {
      if (count(bin) > 0) {
        data(bin) = data(bin) / count(bin);
      }
    }
  }
}

} // namespace proj_details

warper_product_1d_interp_nearest_t
make_warper_from_proj(std::function<dcomplex(std::vector<double>)> const &warped_integrand, int const order,
                      int const num_bins, int const nr_samples, const double sigma) {

  warper_product_1d_interp_nearest_t warper;
  nda::vector<double> time_pts(num_bins + 1);
  nda::vector<double> f1_pts(num_bins + 1);
  std::vector<binner::binner_t<1, 0, double>> xi;
  nda::vector<double> sigmas(order);

  double const delta = 1.0 / num_bins;
  for (int i = 0; i < num_bins + 1; i++) {
    time_pts[i] = i * delta;
  }
  time_pts[num_bins] = 1.;

  for (int axis = 0; axis < order; axis++) {
    xi.emplace_back(binner::binner_t<1, 0, double>({std::make_tuple(0., 1., num_bins)}));
  }

  proj_details::gather_data(warped_integrand, xi, nr_samples, order);

  for (int axis = 0; axis < order; axis++) {
    sigmas[axis] = sigma / xi[axis].get_bin_size(); // change of unit

    auto smoothed_proj = proj_details::kernel_smoothing(xi[axis].get_data(), sigmas[axis]);

    f1_pts[0] = 0;
    for (int i = 0; i < num_bins; i++) {
      f1_pts[i + 1] = 2 * smoothed_proj(i) - f1_pts[i];
    }
    warper.emplace_back(warper_product_1d_simple_interp_nearest_t(time_pts, f1_pts));
  }

  return warper;
}

} // namespace keldy::warpers
