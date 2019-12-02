#include <gtest/gtest.h>
#include <keldy/interfaces/cuba_interface.hpp>
#include <triqs/test_tools/gfs.hpp>
#include <itertools/itertools.hpp>
#include <boost/math/special_functions/erf.hpp>

using namespace keldy;

// Note: These are just tests of the wrapping, not integragtor benchmarks.

double gaussian(std::vector<double> times) {
  double result = 1.0;
  for (auto const &[i, t] : itertools::enumerate(times)) {
    double two_var = 2.0 / std::pow(2, i);
    result *= std::exp(-std::pow(t - 0.5, 2) / two_var) / std::sqrt(M_PI * two_var);
  }
  return result;
}

double gaussian_hypercube_analytic(int dim) {
  double analytic_result = 1.0;
  for (int i = 0; i < dim; i++) {
    analytic_result *= boost::math::erf(std::sqrt(std::pow(2.0, i - 3)));
  }
  return analytic_result;
}

TEST(cuba_wrapper, cuhre_gaussian) { // NOLINT
  cuba_common_param p;
  p.error_eps_rel = 1e-3;
  p.error_eps_abs = 1e-4;

  std::vector<int> dim_v{4, 6, 10};
  for (auto dim : dim_v) {
    auto analytic = gaussian_hypercube_analytic(dim);

    cuba_wrapper cuba{gaussian, dim, p};
    cuba.run_cuhre(0);
    auto o = cuba.get_output();

    // Basic Sanity Checks
    EXPECT_GE(o.n_regions, 1);
    EXPECT_GE(o.n_eval, p.min_number_evaluations);
    EXPECT_GE(o.chi_sq_prob, 0.0);

    EXPECT_TRUE(o.error_flag == 0); // Successful Termination

    double rerr_estimate = p.error_eps_rel * o.result;
    if (rerr_estimate < p.error_eps_abs) {
      EXPECT_LE(std::abs(analytic - o.result), 2 * p.error_eps_abs);
    } else {
      EXPECT_LE(std::abs(analytic - o.result), 2 * rerr_estimate);
    }

    EXPECT_LE(std::abs(analytic - o.result), 2 * o.error); // Good Error Estimate
  }
}

TEST(cuba_wrapper, vegas_gaussian) { // NOLINT
  cuba_common_param p;
  p.error_eps_rel = 1e-3;
  p.error_eps_abs = 1e-4;

  cuba_vegas_param p_v;

  std::vector<int> dim_v{4, 6, 10};
  for (auto dim : dim_v) {
    auto analytic = gaussian_hypercube_analytic(dim);

    cuba_wrapper cuba{gaussian, dim, p};
    cuba.run_vegas(p_v);
    auto o = cuba.get_output();

    // Basic Sanity Checks
    EXPECT_EQ(o.n_regions, 0); // No Regions in Vegas algorithm
    EXPECT_GE(o.n_eval, p.min_number_evaluations);
    EXPECT_GE(o.chi_sq_prob, 0.0);

    EXPECT_TRUE(o.error_flag == 0); // Successful Termination

    double rerr_estimate = p.error_eps_rel * o.result;
    if (rerr_estimate < p.error_eps_abs) {
      EXPECT_LE(std::abs(analytic - o.result), 2 * p.error_eps_abs);
    } else {
      EXPECT_LE(std::abs(analytic - o.result), 2 * rerr_estimate);
    }

    EXPECT_LE(std::abs(analytic - o.result), 2 * o.error); // Good Error Estimate
  }
}

TEST(cuba_wrapper, suave_gaussian) { // NOLINT
  cuba_common_param p;
  p.error_eps_rel = 1e-3;
  p.error_eps_abs = 1e-4;

  cuba_suave_param p_s;
  p_s.n_new_evals_each_subdivision = 5000;
  p_s.n_min_samples_region_threashold = 500;
  p_s.flatness_parameter_p = 50.0;

  std::vector<int> dim_v{4, 6, 10};
  for (auto dim : dim_v) {
    auto analytic = gaussian_hypercube_analytic(dim);

    cuba_wrapper cuba{gaussian, dim, p};
    cuba.run_suave(p_s);
    auto o = cuba.get_output();

    // Basic Sanity Checks
    EXPECT_GE(o.n_regions, 1);
    EXPECT_GE(o.n_eval, p.min_number_evaluations);
    EXPECT_GE(o.chi_sq_prob, 0.0);

    ASSERT_TRUE(o.error_flag == 0); // Successful Termination

    double rerr_estimate = p.error_eps_rel * std::abs(o.result);
    if (p.error_eps_abs > rerr_estimate) {
      EXPECT_LE(std::abs(analytic - o.result), 2 * p.error_eps_abs);
    } else {
      EXPECT_LE(std::abs(analytic - o.result), 2 * rerr_estimate);
    }

    EXPECT_LE(std::abs(analytic - o.result), 2 * o.error); // Good Error Estimate
  }
}

MAKE_MAIN; // NOLINT
