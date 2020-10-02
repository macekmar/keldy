/*******************************************************************************
 *
 * keldy
 *
 * Copyright (c) 2019-2020 The Simons Foundation
 *   authors: Philipp Dumitrescu, Corentin Bertrand
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

#include "model.hpp"
#include "../contour_integral.hpp"
#include "keldy/common.hpp"
#include <limits>
#include <triqs/gfs.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <itertools/omp_chunk.hpp>

namespace keldy::impurity_oneband {

void h5_write(h5::group &h5group, std::string const &subgroup_name, model_param_t const &c) {
  auto grp = h5group.create_group(subgroup_name);
  h5_write(grp, "beta", c.beta);
  h5_write(grp, "bias_V_left", c.bias_V_left);
  h5_write(grp, "bias_V_right", c.bias_V_right);
  h5_write(grp, "eps_d", c.eps_d);
  h5_write(grp, "Gamma", c.Gamma);
  h5_write(grp, "alpha", c.alpha);
  h5_write(grp, "half_bandwidth", c.half_bandwidth);
  h5_write(grp, "bath_type", c.bath_type);
  h5_write(grp, "time_max", c.time_max);
  h5_write(grp, "nr_time_points_gf", c.nr_time_points_gf);
  h5_write(grp, "ft_method", c.ft_method);
}

void h5_read(h5::group &h5group, std::string const &subgroup_name, model_param_t &c) {
  auto grp = h5group.open_group(subgroup_name);
  h5_read(grp, "beta", c.beta);
  h5_read(grp, "bias_V_left", c.bias_V_left);
  h5_read(grp, "bias_V_right", c.bias_V_right);
  h5_read(grp, "eps_d", c.eps_d);
  h5_read(grp, "Gamma", c.Gamma);
  h5_read(grp, "alpha", c.alpha);
  h5_read(grp, "half_bandwidth", c.half_bandwidth);
  h5_read(grp, "bath_type", c.bath_type);
  h5_read(grp, "time_max", c.time_max);
  h5_read(grp, "nr_time_points_gf", c.nr_time_points_gf);
  h5_read(grp, "ft_method", c.ft_method);
}

// *****************************************************************************

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
  using namespace std::complex_literals;
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
    bath_hybrid_R_left_ = [Gamma = param_.Gamma, half_bandwidth = param_.half_bandwidth](dcomplex omega) -> dcomplex {
      if (std::imag(omega) != 0.) {
        TRIQS_RUNTIME_ERROR << "Semicircular hybridization function out of the real axis is not implemented.";
      }
      if (std::abs(std::real(omega)) < half_bandwidth) {
        return 0.5 * (Gamma / half_bandwidth)
           * (omega - 1.0i * std::sqrt(half_bandwidth * half_bandwidth - omega * omega));
      }
      auto sgn_omega = (1 - 2 * int(std::signbit(std::real(omega))));
      return 0.5 * (Gamma / half_bandwidth)
         * (omega - sgn_omega * std::sqrt(omega * omega - half_bandwidth * half_bandwidth));
    };
    bath_hybrid_R_right_ = bath_hybrid_R_left_;

  } else if (param_.bath_type == "flatband") {
    bath_hybrid_R_left_ = [Gamma = param_.Gamma]([[maybe_unused]] dcomplex omega) -> dcomplex {
      return - 1.0i * Gamma / 2.;
    };
    bath_hybrid_R_right_ = bath_hybrid_R_left_;

    /// half_bandwidth is here interpreted as half width at half maximum.
  } else if (param_.bath_type == "lorentzian") {
    bath_hybrid_R_left_ = [Gamma = param_.Gamma, half_bandwidth = param_.half_bandwidth](dcomplex omega) -> dcomplex {
      return Gamma * half_bandwidth / 2. / (omega + 1.0i * half_bandwidth);
    };
    bath_hybrid_R_right_ = bath_hybrid_R_left_;

  } else {
    TRIQS_RUNTIME_ERROR << "bath_type not defined";
  }
}

void h5_write(h5::group &h5group, std::string const &subgroup_name, g0_model_omega const &c) {
  auto grp = h5group.create_group(subgroup_name);
  h5_write(grp, "param_", c.param_);
}

void h5_read(h5::group &h5group, std::string const &subgroup_name, g0_model_omega &c) {
  auto grp = h5group.open_group(subgroup_name);
  model_param_t param_temp;
  h5_read(grp, "param_", param_temp);
  c = g0_model_omega{param_temp};
}

// *****************************************************************************

g0_model::g0_model(g0_model_omega model_omega_, bool make_dot_lead_)
   : model_omega(std::move(model_omega_)), make_dot_lead(make_dot_lead_) {

  auto param_ = model_omega.get_param();

  if (make_dot_lead && (param_.bath_type == "flatband")) {
    TRIQS_RUNTIME_ERROR << "The dot-lead g0 for flatband bath is singular and not implemented accurately.";
  }

  if (param_.ft_method == "fft") {
    make_g0_by_fft();

  } else if (param_.ft_method == "contour") {
    auto half_bias = std::abs(param_.bias_V_left - param_.bias_V_right) / 2;
    if (param_.bath_type == "semicircle") { // semicircular DOS has a bounded support
      if (half_bias < param_.half_bandwidth) {
        make_g0_by_finite_contour({-param_.half_bandwidth, -half_bias, half_bias, param_.half_bandwidth});
      } else {
        make_g0_by_finite_contour({-param_.half_bandwidth, param_.half_bandwidth});
      }
    } else {
      // TODO: use qagp ?
      double margin = std::abs(std::max(param_.Gamma, 1. / param_.beta));
      double left_turn_pt = std::min(-half_bias, param_.eps_d) - margin;
      double right_turn_pt = std::max(half_bias, param_.eps_d) + margin;
      make_g0_by_contour(left_turn_pt, right_turn_pt);
    }

  } else if (param_.ft_method == "analytic") {
    if (param_.bath_type == "flatband") {
      make_flat_band_analytic();
    } else {
      TRIQS_RUNTIME_ERROR << "analytic ft only works for flatband";
    }

  } else {
    TRIQS_RUNTIME_ERROR << "ft_method '" << param_.ft_method << "' is not defined";
  }
}

void g0_model::make_g0_by_fft() {
  using namespace std::complex_literals;
  using namespace triqs::gfs;

  auto param_ = model_omega.get_param();

  auto time_mesh = gf_mesh<retime>({-param_.time_max, param_.time_max, param_.nr_time_points_gf});
  auto freq_mesh = make_adjoint_mesh(time_mesh);

  gf<refreq, matrix_valued> g0_lesser_omega{freq_mesh, {2, 2}};
  gf<refreq, matrix_valued> g0_greater_omega{freq_mesh, {2, 2}};

  for (auto omega : freq_mesh) {
    dcomplex omegaj = omega + 0.0i;

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
 * Compute the Fourier transform of g^{</>}(\omega) by integrating along a finite interval.
 *
 * This method is suitable only for bath dos with a finite support, like semi-circular.
 * omega is real in this method.
 */
void g0_model::make_g0_by_finite_contour(std::vector<double> pts) {
  using namespace triqs::gfs;
  using namespace boost::math::double_constants;
  using namespace std::complex_literals;

  // sort and remove duplicates from pts
  std::sort(pts.begin(), pts.end());
  pts.erase(std::unique(pts.begin(), pts.end()), pts.end());

  auto param_ = model_omega.get_param();

  auto time_mesh = gf_mesh<retime>({-param_.time_max, param_.time_max, param_.nr_time_points_gf});

  gf<retime, matrix_valued> g0_lesser_time{time_mesh, {2, 2}};
  gf<retime, matrix_valued> g0_greater_time{time_mesh, {2, 2}};
  lesser_ft_error = {time_mesh, {2, 2}};
  greater_ft_error = {time_mesh, {2, 2}};

  auto gsl_default_handler = gsl_set_error_handler_off();
  double const abstol = 1e-14;
  double const reltol = 1e-8;

#pragma omp parallel
  {
    details::gsl_integration_wrapper worker{1000};
    dcomplex result;
    dcomplex abserr;

    for (const auto t : itertools::omp_chunk(time_mesh)) {

      // g_lesser:
      auto integrand_dot_lesser = [t, model = this->model_omega](double omega) -> dcomplex {
        return std::exp(-1.0i * omega * t) * model.g0_dot_lesser(omega) / (2 * pi);
      };

      std::tie(result, abserr) = worker.qagp(integrand_dot_lesser, pts, abstol, reltol);
      g0_lesser_time[t](0, 0) = result;
      lesser_ft_error[t](0, 0) = abserr;

      if (make_dot_lead) {
        auto integrand_lead_lesser = [t, model = this->model_omega](double omega) -> dcomplex {
          return std::exp(-1.0i * omega * t) * model.g0_rightlead_dot_lesser(omega) / (2 * pi);
        };

        std::tie(result, abserr) = worker.qagp(integrand_lead_lesser, pts, abstol, reltol);
        g0_lesser_time[t](1, 0) = result;
        lesser_ft_error[t](1, 0) = abserr;
      }

      // g_greater
      auto integrand_dot_greater = [t, model = this->model_omega](double omega) -> dcomplex {
        return std::exp(-1.0i * omega * t) * model.g0_dot_greater(omega) / (2 * pi);
      };

      std::tie(result, abserr) = worker.qagp(integrand_dot_greater, pts, abstol, reltol);
      g0_greater_time[t](0, 0) = result;
      greater_ft_error[t](0, 0) = abserr;

      if (make_dot_lead) {
        auto integrand_lead_greater = [t, model = this->model_omega](double omega) -> dcomplex {
          return std::exp(-1.0i * omega * t) * model.g0_rightlead_dot_greater(omega) / (2 * pi);
        };

        std::tie(result, abserr) = worker.qagp(integrand_lead_greater, pts, abstol, reltol);
        g0_greater_time[t](1, 0) = result;
        greater_ft_error[t](1, 0) = abserr;
      }
    }

    // TODO: make use of t <-> -t symmetry to reduce calculations
  }

  // GSL error estimations
  //lesser_ft_error(0, 1) = lesser_ft_error(1, 0);
  //greater_ft_error(0, 1) = greater_ft_error(1, 0);

  // reset default GSL error handler
  gsl_set_error_handler(gsl_default_handler);

  // Spin up and down are currently identical
  g0_lesser = make_block_gf<retime, matrix_valued>({"up", "down"}, {g0_lesser_time, g0_lesser_time});
  g0_greater = make_block_gf<retime, matrix_valued>({"up", "down"}, {g0_greater_time, g0_greater_time});
}

/*
 * Compute the Fourier transform of g^{</>}(\omega) by integrating along a suitable contour in the complex plane.
 *
 * Use std::imag and std::conj carefully, omega is complex here!!
 */
void g0_model::make_g0_by_contour(double left_turn_pt, double right_turn_pt) {
  using namespace triqs::gfs;
  using namespace boost::math::double_constants;
  using namespace std::complex_literals;

  auto param_ = model_omega.get_param();

  if (!(left_turn_pt < -std::abs(param_.bias_V_left) && right_turn_pt > std::abs(param_.bias_V_right))) {
    TRIQS_RUNTIME_ERROR << "Contour is wrong regarding Fermi functions poles.";
  }

  auto time_mesh = gf_mesh<retime>({-param_.time_max, param_.time_max, param_.nr_time_points_gf});

  gf<retime, matrix_valued> g0_lesser_time{time_mesh, {2, 2}};
  gf<retime, matrix_valued> g0_greater_time{time_mesh, {2, 2}};
  lesser_ft_error = {time_mesh, {2, 2}};
  greater_ft_error = {time_mesh, {2, 2}};

  auto gsl_default_handler = gsl_set_error_handler_off();

#pragma omp parallel
  {
    contour_integration_t worker(left_turn_pt, right_turn_pt);

    for (const auto t : itertools::omp_chunk(time_mesh)) {
      int sign_of_t = (t > 0) - (t < 0);
      dcomplex direc_1 = (t == 0) ? -1 : -1.0i * sign_of_t;
      dcomplex direc_2 = (param_.beta == std::numeric_limits<double>::infinity()) ? 1. : param_.beta - t * 1.0i;
      direc_2 /= std::abs(direc_2);

      // g_lesser:
      auto integrand_dot_lesser = [t, model = this->model_omega](dcomplex omega) -> dcomplex {
        return std::exp(-1.0i * omega * t) * model.g0_dot_lesser(omega) / (2 * pi);
      };

      worker.integrate(integrand_dot_lesser, direc_1, direc_2);
      g0_lesser_time[t](0, 0) = worker.get_result();
      lesser_ft_error[t](0, 0) = worker.get_abserr();

      if (make_dot_lead) {
        auto integrand_lead_lesser = [t, model = this->model_omega](dcomplex omega) -> dcomplex {
          return std::exp(-1.0i * omega * t) * model.g0_rightlead_dot_lesser(omega) / (2 * pi);
        };

        worker.integrate(integrand_lead_lesser, direc_1, direc_2);
        g0_lesser_time[t](1, 0) = worker.get_result();
        lesser_ft_error[t](1, 0) = worker.get_abserr();
      }

      // g_greater
      auto integrand_dot_greater = [t, model = this->model_omega](dcomplex omega) -> dcomplex {
        return std::exp(-1.0i * omega * t) * model.g0_dot_greater(omega) / (2 * pi);
      };

      worker.integrate(integrand_dot_greater, -std::conj(direc_2), -std::conj(direc_1));
      g0_greater_time[t](0, 0) = worker.get_result();
      greater_ft_error[t](0, 0) = worker.get_abserr();

      if (make_dot_lead) {
        auto integrand_lead_greater = [t, model = this->model_omega](dcomplex omega) -> dcomplex {
          return std::exp(-1.0i * omega * t) * model.g0_rightlead_dot_greater(omega) / (2 * pi);
        };

        worker.integrate(integrand_lead_greater, -std::conj(direc_2), -std::conj(direc_1));
        g0_greater_time[t](1, 0) = worker.get_result();
        greater_ft_error[t](1, 0) = worker.get_abserr();
      }
    }

    // TODO: make use of t <-> -t symmetry to reduce calculations
  }

  // GSL error estimations
  //lesser_ft_error(0, 1) = lesser_ft_error(1, 0);
  //greater_ft_error(0, 1) = greater_ft_error(1, 0);

  // reset default GSL error handler
  gsl_set_error_handler(gsl_default_handler);

  // Spin up and down are currently identical
  g0_lesser = make_block_gf<retime, matrix_valued>({"up", "down"}, {g0_lesser_time, g0_lesser_time});
  g0_greater = make_block_gf<retime, matrix_valued>({"up", "down"}, {g0_greater_time, g0_greater_time});
}

void g0_model::make_flat_band_analytic() {
  using namespace std::complex_literals;

  auto param_ = model_omega.get_param();
  if (param_.beta != std::numeric_limits<double>::infinity() || param_.bias_V_left != 0. || param_.bias_V_right != 0.0 || param_.eps_d != 0.) {
    TRIQS_RUNTIME_ERROR
       << "Analytic flatband covers only the following parameters: beta=infinity (zero temperature), bias_V=0, eps_d=0.";
  }
  if (make_dot_lead) {
    TRIQS_RUNTIME_ERROR << "No analytic formula for dot-lead Green's functions.";
  }

  auto const time_mesh = gf_mesh<retime>({-param_.time_max, param_.time_max, param_.nr_time_points_gf});
  gf<retime, matrix_valued> g0_lesser_up{time_mesh, {2, 2}};
  gf<retime, matrix_valued> g0_greater_up{time_mesh, {2, 2}};

  auto const g0_lesser_values = [Gamma = param_.Gamma](time_real_t time) -> dcomplex {
    using namespace boost::math;
    using namespace boost::math::double_constants;

    if (time == 0.0) {
      return 0.5i;
    }
    auto Gt = Gamma * time;
    auto real_part = (std::exp(Gt) * boost::math::expint(-Gt) - std::exp(-Gt) * boost::math::expint(Gt)) / (2 * pi);
    return real_part + 0.5i * std::exp(-std::abs(Gt));
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

using namespace triqs::arrays;

g0_model::g0_model(model_param_t const &parameters, int n, array<dcomplex, 3> const &g0_lesser_data, array<dcomplex, 3> const &g0_greater_data) {
  // Creates g0(t) from data
  auto param_ = parameters;
  auto const time_mesh = gf_mesh<retime>({-param_.time_max, param_.time_max, param_.nr_time_points_gf});

  gf<retime, matrix_valued> g0_lesser_up{time_mesh, {n, n}};
  gf<retime, matrix_valued> g0_greater_up{time_mesh, {n, n}};

  // Loop to copy data
  for (auto const [i, t] : itertools::enumerate(time_mesh)) {
    g0_lesser_up[t] = g0_lesser_data(i,range(),range());
    g0_greater_up[t] = g0_greater_data(i,range(),range());
  }
  // Since Spin up and down are currently identical
  g0_lesser = make_block_gf<retime, matrix_valued>({"up", "down"}, {g0_lesser_up, g0_lesser_up});
  g0_greater = make_block_gf<retime, matrix_valued>({"up", "down"}, {g0_greater_up, g0_greater_up});
}

g0_model::g0_model(g0_model_omega model_omega_, 
                   gf<retime, matrix_valued> g0_lesser_up, 
                   gf<retime, matrix_valued> g0_greater_up) : model_omega(std::move(model_omega_)), make_dot_lead(true) {
  // Since Spin up and down are currently identical
  g0_lesser = make_block_gf<retime, matrix_valued>({"up", "down"}, {g0_lesser_up, g0_lesser_up});
  g0_greater = make_block_gf<retime, matrix_valued>({"up", "down"}, {g0_greater_up, g0_greater_up});
}

g0_model::g0_model(model_param_t const &parameters, gf_mesh<refreq> freq_mesh, 
                    gf<refreq, matrix_valued> g0_lesser_omega, 
                    gf<refreq, matrix_valued> g0_greater_omega) : model_omega(parameters) {

  // model_omega = g0_model_omega(parameters);              

  auto time_mesh = make_adjoint_mesh(freq_mesh);

  gf<retime, matrix_valued> g0_lesser_up = make_gf_from_fourier(g0_lesser_omega, time_mesh);
  gf<retime, matrix_valued> g0_greater_up = make_gf_from_fourier(g0_greater_omega, time_mesh);

  // Since Spin up and down are currently identical
  g0_lesser = make_block_gf<retime, matrix_valued>({"up", "down"}, {g0_lesser_up, g0_lesser_up});
  g0_greater = make_block_gf<retime, matrix_valued>({"up", "down"}, {g0_greater_up, g0_greater_up});
}


} // namespace keldy::impurity_oneband
