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

#include "model.hpp"
#include "../contour_integral.hpp"
#include "keldy/common.hpp"
#include <limits>
#include <triqs/gfs.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <itertools/omp_chunk.hpp>

namespace keldy::impurity_oneband {

bool operator<(gf_index_t const &a, gf_index_t const &b) {
  // First: order on the time contour
  int contour_3way = compare_3way(a.contour, b.contour);
  if (contour_3way != 0) {
    return contour_3way < 0;
  }
  // Second: order orbital indices
  if (a.orbital != b.orbital) {
    return (a.orbital < b.orbital);
  }
  // Third: order spin
  return a.spin < b.spin;
}

// *****************************************************************************

g0_model_omega::g0_model_omega(model_param_t const &parameters) : param_(parameters) {
  if (param_.beta < 0) {
    std::cout << "WARNING: beta < 0. Interpret this as beta = +infinity." << std::endl;
    param_.beta = std::numeric_limits<double>::infinity();
  }

  // Ensure t=0 is in the mesh
  if (param_.nr_time_points_gf % 2 == 0) {
    std::cout << "WARNING: number of time points has been increased to make it odd." << std::endl;
    param_.nr_time_points_gf += 1;
  }

  if (param_.bath_type == "semicircle") {
    bath_hybrid_R_left_ = [&Gamma = param_.Gamma](dcomplex omega) -> dcomplex {
      auto omega2 = omega / 2.;
      if (std::abs(std::real(omega2)) < 1) {
        return (Gamma / 2) * (omega2 - 1_j * std::sqrt(1 - omega2 * omega2));
      }
      auto sgn_omega = (1 - 2 * int(std::signbit(std::real(omega2))));
      return (Gamma / 2) * (omega2 - sgn_omega * std::sqrt(omega2 * omega2 - 1));
    };
    bath_hybrid_R_right_ = bath_hybrid_R_left_;

  } else if (param_.bath_type == "flatband") {
    bath_hybrid_R_left_ = [&Gamma = param_.Gamma]([[maybe_unused]] dcomplex omega) -> dcomplex {
      return -1_j * Gamma / 2.;
    };
    bath_hybrid_R_right_ = bath_hybrid_R_left_;

  } else {
    TRIQS_RUNTIME_ERROR << "bath_type not defined";
  }
}

// *****************************************************************************

g0_model::g0_model(g0_model_omega model_omega_, bool make_dot_lead_)
   : model_omega(std::move(model_omega_)), make_dot_lead(make_dot_lead_) {

  if (model_omega.param_.ft_method == "fft") {
    make_g0_by_fft();

  } else if (model_omega.param_.ft_method == "contour") {
    double margin = std::abs(std::max(model_omega.param_.Gamma, 1. / model_omega.param_.beta));
    double left_turn_pt = std::min(-std::abs(model_omega.param_.bias_V / 2), model_omega.param_.eps_d) - margin;
    double right_turn_pt = std::max(std::abs(model_omega.param_.bias_V / 2), model_omega.param_.eps_d) + margin;
    make_g0_by_contour(left_turn_pt, right_turn_pt);

  } else if (model_omega.param_.ft_method == "analytic") {
    if (model_omega.param_.bath_type == "flatband") {
      make_flat_band_analytic();
    } else {
      TRIQS_RUNTIME_ERROR << "analytic ft only works for flatband";
    }

  } else {
    TRIQS_RUNTIME_ERROR << "ft_method not defined";
  }
}

void g0_model::make_g0_by_fft() {
  using namespace triqs::gfs;
  auto time_mesh =
     gf_mesh<retime>({-model_omega.param_.time_max, model_omega.param_.time_max, model_omega.param_.nr_time_points_gf});
  auto freq_mesh = make_adjoint_mesh(time_mesh);

  gf<refreq, matrix_valued> g0_lesser_omega{freq_mesh, {2, 2}};
  gf<refreq, matrix_valued> g0_greater_omega{freq_mesh, {2, 2}};

  for (auto omega : freq_mesh) {
    dcomplex omegaj = omega + 0_j;

    g0_lesser_omega[omega](0, 0) = model_omega.g0_dot_lesser(omegaj);
    g0_greater_omega[omega](0, 0) = model_omega.g0_dot_greater(omegaj);

    if (make_dot_lead) {
      /// right lead only
      g0_lesser_omega[omega](1, 0) = model_omega.g0_rightlead_dot_lesser(omegaj);
      g0_greater_omega[omega](1, 0) = model_omega.g0_rightlead_dot_greater(omegaj);
      g0_lesser_omega[omega](0, 1) = -std::conj(g0_lesser_omega[omega](1, 0));
      g0_greater_omega[omega](0, 1) = -std::conj(g0_greater_omega[omega](1, 0));
    }
  }

  gf<retime, matrix_valued> g0_lesser_up = make_gf_from_fourier(g0_lesser_omega, time_mesh);
  gf<retime, matrix_valued> g0_greater_up = make_gf_from_fourier(g0_greater_omega, time_mesh);

  // Since Spin up and down are currently identical
  g0_lesser = make_block_gf<retime, matrix_valued>({"up", "down"}, {g0_lesser_up, g0_lesser_up});
  g0_greater = make_block_gf<retime, matrix_valued>({"up", "down"}, {g0_greater_up, g0_greater_up});
}

/*
 * Compute the Fourier transform of g^{</>}(\omega) by integrating along a suitable contour in the complex plane.
 *
 * Use std::imag and std::conj carefully, omega is complex here!!
 */
void g0_model::make_g0_by_contour(double left_turn_pt, double right_turn_pt) {
  using namespace triqs::gfs;
  using namespace boost::math::double_constants;

  if (!(left_turn_pt < -std::abs(model_omega.param_.bias_V / 2)
        && right_turn_pt > std::abs(model_omega.param_.bias_V / 2))) {
    TRIQS_RUNTIME_ERROR << "Contour is wrong regarding Fermi functions poles.";
  }

  auto time_mesh =
     gf_mesh<retime>({-model_omega.param_.time_max, model_omega.param_.time_max, model_omega.param_.nr_time_points_gf});

  gf<retime, matrix_valued> g0_lesser_time{time_mesh, {2, 2}};
  gf<retime, matrix_valued> g0_greater_time{time_mesh, {2, 2}};

  auto gsl_default_handler = gsl_set_error_handler_off();

#pragma omp parallel
  {
    contour_integration_t worker(left_turn_pt, right_turn_pt);

    for (const auto t : itertools::omp_chunk(time_mesh)) {
      int sign_of_t = (t > 0) - (t < 0);
      dcomplex direc_1 = (t == 0) ? -1 : -1_j * sign_of_t;
      dcomplex direc_2 =
         (model_omega.param_.beta == std::numeric_limits<double>::infinity()) ? 1. : model_omega.param_.beta - t * 1_j;
      direc_2 /= std::abs(direc_2);

      // g_lesser:
      auto integrand_dot_lesser = [t, model = this->model_omega](dcomplex omega) -> dcomplex {
        return std::exp(-1_j * omega * t) * model.g0_dot_lesser(omega) / (2 * pi);
      };

      worker.integrate(integrand_dot_lesser, direc_1, direc_2);
      g0_lesser_time[t](0, 0) = worker.get_result();
      lesser_ft_error(0, 0) += worker.get_abserr_sqr();

      if (make_dot_lead) {
        auto integrand_lead_lesser = [t, model = this->model_omega](dcomplex omega) -> dcomplex {
          return std::exp(-1_j * omega * t) * model.g0_rightlead_dot_lesser(omega) / (2 * pi);
        };

        worker.integrate(integrand_lead_lesser, direc_1, direc_2);
        g0_lesser_time[t](1, 0) = worker.get_result();
        lesser_ft_error(1, 0) += worker.get_abserr_sqr();
      }

      // g_greater
      auto integrand_dot_greater = [t, model = this->model_omega](dcomplex omega) -> dcomplex {
        return std::exp(-1_j * omega * t) * model.g0_dot_greater(omega) / (2 * pi);
      };

      worker.integrate(integrand_dot_greater, -std::conj(direc_2), -std::conj(direc_1));
      g0_greater_time[t](0, 0) = worker.get_result();
      greater_ft_error(0, 0) += worker.get_abserr_sqr();

      if (make_dot_lead) {
        auto integrand_lead_greater = [t, model = this->model_omega](dcomplex omega) -> dcomplex {
          return std::exp(-1_j * omega * t) * model.g0_rightlead_dot_greater(omega) / (2 * pi);
        };

        worker.integrate(integrand_lead_greater, -std::conj(direc_2), -std::conj(direc_1));
        g0_greater_time[t](1, 0) = worker.get_result();
        greater_ft_error(1, 0) += worker.get_abserr_sqr();
      }
    }

    // TODO: make use of t <-> -t symmetry to reduce calculations
  }

  // GSL error estimations
  //lesser_ft_error(0, 1) = lesser_ft_error(1, 0);
  lesser_ft_error = sqrt(lesser_ft_error / time_mesh.size());
  //greater_ft_error(0, 1) = greater_ft_error(1, 0);
  greater_ft_error = sqrt(greater_ft_error / time_mesh.size());

  // reset default GSL error handler
  gsl_set_error_handler(gsl_default_handler);

  // Spin up and down are currently identical
  g0_lesser = make_block_gf<retime, matrix_valued>({"up", "down"}, {g0_lesser_time, g0_lesser_time});
  g0_greater = make_block_gf<retime, matrix_valued>({"up", "down"}, {g0_greater_time, g0_greater_time});
}

void g0_model::make_flat_band_analytic() {
  if (model_omega.param_.beta != std::numeric_limits<double>::infinity() || model_omega.param_.bias_V != 0. || model_omega.param_.eps_d != 0.) {
    TRIQS_RUNTIME_ERROR
       << "Analytic flatband covers only the following parameters: beta=infinity (zero temperature), bias_V=0, eps_d=0.";
  }
  if (make_dot_lead) {
    TRIQS_RUNTIME_ERROR << "No analytic formula for dot-lead Green's functions.";
  }

  auto const time_mesh = gf_mesh<retime>({-model_omega.param_.time_max, model_omega.param_.time_max, model_omega.param_.nr_time_points_gf});
  gf<retime, matrix_valued> g0_lesser_up{time_mesh, {2, 2}};
  gf<retime, matrix_valued> g0_greater_up{time_mesh, {2, 2}};

  auto const g0_lesser_values = [Gamma = model_omega.param_.Gamma](time_real_t time) -> dcomplex {
    using namespace boost::math;
    using namespace boost::math::double_constants;

    if (time == 0.0) {
      return 0.5_j;
    }
    auto Gt = Gamma * time;
    auto real_part = (std::exp(Gt) * expint(-Gt) - std::exp(-Gt) * expint(Gt)) / (2 * pi);
    return real_part + 0.5_j * std::exp(-std::abs(Gt));
  };

  dcomplex val = 0.;
  for (auto t : time_mesh) {
    val = g0_lesser_values(t);
    g0_lesser_up[t](0, 0) = val;
    g0_greater_up[t](0, 0) = std::conj(val);
  }

  // Since Spin up and down are currently identical
  g0_lesser = make_block_gf<retime, matrix_valued>({"up", "down"}, {g0_lesser_up, g0_lesser_up});
  g0_greater = make_block_gf<retime, matrix_valued>({"up", "down"}, {g0_greater_up, g0_greater_up});
}

} // namespace keldy::impurity_oneband
