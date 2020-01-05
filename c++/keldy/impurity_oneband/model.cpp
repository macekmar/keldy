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
#include <limits>
#include <triqs/gfs.hpp>
#include <boost/math/special_functions/expint.hpp>

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

// *****

g0_model::g0_model(model_param_t const &parameters, bool with_leads) : param_(parameters), contain_leads(with_leads) {
  if (param_.beta < 0) {
    std::cout << "WARNING: beta < 0. Interpret this as beta = +infinity." << std::endl;
    param_.beta = std::numeric_limits<double>::infinity();
  }

  if (param_.bath_type == "semicircle_fft") {
    bath_hybrid_R_left = [Gamma = param_.Gamma](dcomplex omega) -> dcomplex {
      auto omega2 = omega / 2.;
      if (std::abs(std::real(omega2)) < 1) {
        return (Gamma / 2) * (omega2 - 1_j * std::sqrt(1 - omega2 * omega2));
      }
      auto sgn_omega = (1 - 2 * int(std::signbit(std::real(omega2))));
      return (Gamma / 2) * (omega2 - sgn_omega * std::sqrt(omega2 * omega2 - 1));
    };
    bath_hybrid_R_right_ = bath_hybrid_R_left_;
    make_g0_by_fft();

  } else if (param_.bath_type == "flatband_fft") {
    bath_hybrid_R_left_ = [Gamma = param_.Gamma]([[maybe_unused]] dcomplex omega) { return -1_j * Gamma / 2.; };
    bath_hybrid_R_right_ = bath_hybrid_R_left_;
    make_g0_by_fft();

  } else if (param_.bath_type == "flatband_contour") {
    bath_hybrid_R_left_ = [Gamma = param_.Gamma]([[maybe_unused]] dcomplex omega) { return -1_j * Gamma / 2.; };
    bath_hybrid_R_right_ = bath_hybrid_R_left_;

    double margin = std::abs(std::max(param_.Gamma, 1. / param_.beta));
    double left_turn_pt = std::min(-std::abs(param_.bias_V / 2), param_.eps_d) - margin;
    double right_turn_pt = std::max(std::abs(param_.bias_V / 2), param_.eps_d) + margin;

    make_g0_by_contour(left_turn_pt, right_turn_pt);

    //FIXME: the lead-dot green's functions at t=0 is infinite (logarithmic singularity) for flatband. How to treat it correctly?

  } else if (param_.bath_type == "flatband_analytic") {
    // bath_hybrid_R not used to get g(t)
    bath_hybrid_R_left_ = [Gamma = param_.Gamma]([[maybe_unused]] dcomplex omega) { return -1_j * Gamma / 2.; };
    bath_hybrid_R_right_ = bath_hybrid_R_left_;
    make_flat_band_analytic();
  } else {
    TRIQS_RUNTIME_ERROR << "bath_type not defined";
  }
}

void g0_model::make_g0_by_fft() {
  using namespace triqs::gfs;

  auto time_mesh = gf_mesh<retime>({-param_.time_max, param_.time_max, param_.nr_time_points_gf});
  auto freq_mesh = make_adjoint_mesh(time_mesh);

  gf<refreq, matrix_valued> g0_lesser_omega{freq_mesh, {2, 2}};
  gf<refreq, matrix_valued> g0_greater_omega{freq_mesh, {2, 2}};

  for (auto omega : freq_mesh) {
    dcomplex omegaj = omega + 0_j;
    auto R = g0_R_dot(omegaj);
    auto A = g0_A_dot(omegaj);
    auto K = R * (bath_hybrid_K_left(omegaj) + bath_hybrid_K_right(omegaj)) * A;

    g0_lesser_omega[omega](0, 0) = (K - R + A) / 2;
    g0_greater_omega[omega](0, 0) = (K + R - A) / 2;

    if (contain_leads) {
      /// right lead only
      K = bath_hybrid_K_right(omegaj) * A + bath_hybrid_R_right(omegaj) * K;
      R = bath_hybrid_R_right(omegaj) * g0_R_dot(omegaj);
      A = std::conj(R);
      g0_lesser_omega[omega](1, 0) = (K - R + A) / 2;
      g0_greater_omega[omega](1, 0) = (K + R - A) / 2;
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

  if (not(left_turn_pt < -std::abs(param_.bias_V / 2) && right_turn_pt > std::abs(param_.bias_V / 2))) {
    TRIQS_RUNTIME_ERROR << "Contour is wrong regarding Fermi functions poles.";
  }

  auto time_mesh = gf_mesh<retime>({-param_.time_max, param_.time_max, param_.nr_time_points_gf});

  gf<retime, matrix_valued> g0_lesser_time{time_mesh, {2, 2}};
  gf<retime, matrix_valued> g0_greater_time{time_mesh, {2, 2}};

  double const mu_left = -param_.bias_V / 2;
  double const mu_right = +param_.bias_V / 2;

  // auto bath_hybrid_left_K = [this](dcomplex omega) {
  //   return -(2 * n_fermi(omega + param_.bias_V / 2, param_.beta) - 1.)
  //      * (bath_hybrid_R_left(omega) - bath_hybrid_A_left(omega));
  // };

  // auto bath_hybrid_right_K = [this](dcomplex omega) {
  //   return -(2 * n_fermi(omega - param_.bias_V / 2, param_.beta) - 1.)
  //      * (bath_hybrid_R_right(omega) - bath_hybrid_A_right(omega));
  // };


  auto g0_omega = [&](dcomplex omega, bool is_lesser) {
    auto bath_hybrid_left = bath_hybrid_R_left(omega) - bath_hybrid_A_left(omega);
    auto bath_hybrid_right = bath_hybrid_R_right(omega) - bath_hybrid_A_right(omega);
    int sign = (is_lesser) ? +1 : -1;
    return -sign * g0_R_dot(omega)
       * (n_fermi(sign * (omega - mu_left), param_.beta) * bath_hybrid_left
          + n_fermi(sign * (omega - mu_right), param_.beta) * bath_hybrid_right)
       * g0_A_dot(omega);
  };

  auto g0_rightlead_dot_omega = [&](dcomplex omega, bool is_lesser) {
    auto bath_hybrid_right = bath_hybrid_R_right(omega) - bath_hybrid_A_right(omega);
    int sign = (is_lesser) ? +1 : -1;
    return -sign * n_fermi(sign * (omega - mu_right), param_.beta) * bath_hybrid_right * g0_A_dot(omega)
       + bath_hybrid_R_right(omega) * g0_omega(omega, is_lesser);
  };

  contour_integration_t worker(left_turn_pt, right_turn_pt);
  auto gsl_default_handler = gsl_set_error_handler_off();

  for (const auto t : time_mesh) {
    int sign_of_t = (t > 0) - (t < 0);
    dcomplex direc_1 = (t == 0) ? -1 : -1_j * sign_of_t;
    dcomplex direc_2 = (param_.beta == std::numeric_limits<double>::infinity()) ? 1. : param_.beta - t * 1_j;
    direc_2 /= std::abs(direc_2);

    for (const bool is_lesser : {true, false}) {

      auto integrand = [t, &g0_omega, is_lesser](dcomplex omega) {
        return std::exp(-1_j * omega * t) * g0_omega(omega, is_lesser) / (2 * pi);
      };

      auto integrand_lead = [t, &g0_rightlead_dot_omega, is_lesser](dcomplex omega) {
        return std::exp(-1_j * omega * t) * g0_rightlead_dot_omega(omega, is_lesser) / (2 * pi);
      };

      if (is_lesser) {
        worker.integrate(integrand, direc_1, direc_2);
        g0_lesser_time[t](0, 0) = worker.get_result();
        lesser_ft_error(0, 0) += worker.get_abserr_sqr();

        if (contain_leads) {
          worker.integrate(integrand_lead, direc_1, direc_2);
          g0_lesser_time[t](1, 0) = worker.get_result();
          lesser_ft_error(1, 0) += worker.get_abserr_sqr();
        }
      } else {
        worker.integrate(integrand, -std::conj(direc_2), -std::conj(direc_1));
        g0_greater_time[t](0, 0) = worker.get_result();
        greater_ft_error(0, 0) += worker.get_abserr_sqr();

        if (contain_leads) {
          worker.integrate(integrand_lead, -std::conj(direc_2), -std::conj(direc_1));
          g0_greater_time[t](1, 0) = worker.get_result();
          greater_ft_error(1, 0) += worker.get_abserr_sqr();
        }
      }

      // TODO: make use of t <-> -t symmetry to reduce calculations
    }
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
  if (param_.beta != std::numeric_limits<double>::infinity() || param_.bias_V != 0. || param_.eps_d != 0.) {
    TRIQS_RUNTIME_ERROR
       << "Analytic flatband covers only the following parameters: beta=infinity (zero temperature), bias_V=0, eps_d=0.";
  }
  if (contain_leads) {
    TRIQS_RUNTIME_ERROR << "No analytic formula for dot-lead Green's functions.";
  }

  auto const time_mesh = gf_mesh<retime>({-param_.time_max, param_.time_max, param_.nr_time_points_gf});
  gf<retime, matrix_valued> g0_lesser_up{time_mesh, {2, 2}};
  gf<retime, matrix_valued> g0_greater_up{time_mesh, {2, 2}};

  auto const g0_lesser_values = [Gamma = param_.Gamma](time_real_t time) -> dcomplex {
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
