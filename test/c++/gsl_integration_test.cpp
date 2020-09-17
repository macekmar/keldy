//#include <gtest/gtest.h>
#include <gtest/gtest.h>
#include <keldy/interfaces/gsl_integration_wrap.hpp>
#include <limits>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy::details;
using namespace std::complex_literals;

TEST(gsl_integration, real_qag_lambdas) { // NOLINT
  auto gsl_default_handler = gsl_set_error_handler_off();

  gsl_integration_wrapper worker(100);

  auto function1 = [](double t) -> double { return t * t; };
  auto [result_1, abserr_1] = worker.qag(function1, 1, 2, 1e-10, 1e-7, GSL_INTEG_GAUSS31);
  EXPECT_NEAR(result_1, 2. + 1. / 3., result_1 * 1e-7);

  auto function2 = [](double t) -> double { return 1. / (t * t); };
  auto [result_2, abserr_2] = worker.qag(function2, 1, 2, 1e-10, 1e-7, GSL_INTEG_GAUSS31);
  EXPECT_NEAR(result_2, 1. / 2., result_2 * 1e-7);

  gsl_set_error_handler(gsl_default_handler);
}

TEST(gsl_integration, real_qag_stdfunc) { // NOLINT
  auto gsl_default_handler = gsl_set_error_handler_off();

  gsl_integration_wrapper worker(100);

  std::function<double(double)> function1 = [](double t) -> double { return t * t; };
  auto [result_1, abserr_1] = worker.qag(function1, 1, 2, 1e-10, 1e-7, GSL_INTEG_GAUSS31);
  EXPECT_NEAR(result_1, 2. + 1. / 3., result_1 * 1e-7);

  gsl_set_error_handler(gsl_default_handler);
}

// *********

TEST(gsl_integration, real_qagp_lambdas) { // NOLINT
  auto gsl_default_handler = gsl_set_error_handler_off();

  gsl_integration_wrapper worker(100);

  std::cout << "Function 1" << std::endl;
  auto function1 = [](double t) -> double { return std::sqrt(std::abs(t)); };
  auto [result_1, abserr_1] = worker.qagp(function1, {-1., 0., 1.}, 1e-10, 1e-7);
  EXPECT_NEAR(result_1, 1. + 1. / 3., 1e-7);

  std::cout << "Function 2" << std::endl;
  auto function2 = [](double t) -> double { return 1. / std::sqrt(std::abs(t)); };
  auto [result_2, abserr_2] = worker.qagp(function2, {-1., 0., 1.}, 1e-10, 1e-7);
  EXPECT_NEAR(result_2, 4., 1e-7);

  std::cout << "Function 3" << std::endl;
  auto function3 = [](double t) -> double { return (t < 0.) ? -1. / std::sqrt(-t) : 1. / std::sqrt(t); };
  auto [result_3, abserr_3] = worker.qagp(function3, {-1., 0., 2.}, 1e-10, 1e-7);
  EXPECT_NEAR(result_3, 2. * (std::sqrt(2.) - 1.), 1e-7);

  std::cout << "Function 4" << std::endl;
  auto function4 = [](double t) -> double { return t < 0. ? std::sqrt(1. + t) : 0.5 * std::sqrt(1. - t); };
  auto [result_4, abserr_4] = worker.qagp(function4, {-1., 0., 1.}, 1e-10, 1e-7);
  EXPECT_NEAR(result_4, 1., 1e-7);

  std::cout << "Function 5" << std::endl;
  auto function5 = [](double t) -> double { return t < 0. ? std::sqrt(1. + t) : 0.; };
  auto [result_5, abserr_5] = worker.qagp(function5, {-1., 0., 1.}, 1e-10, 1e-7);
  EXPECT_NEAR(result_5, 2. / 3., 1e-7);

  gsl_set_error_handler(gsl_default_handler);
}

// *********

TEST(gsl_integration, real_qag_si_lambdas) { // NOLINT
  auto gsl_default_handler = gsl_set_error_handler_off();

  gsl_integration_wrapper worker(100);

  double inf = std::numeric_limits<double>::infinity();

  auto function1 = [](double t) -> double { return std::exp(-std::abs(t) * (1 + int(t > 0))); };
  auto [result_1, abserr_1] = worker.qag_si(function1, 0, inf, 1e-10, 1e-7);
  EXPECT_NEAR(result_1, 0.5, 0.5 * 1e-7);
  EXPECT_TRUE(abserr_1 <= 0.5 * 1e-7);

  auto [result_2, abserr_2] = worker.qag_si(function1, -inf, 0, 1e-10, 1e-7);
  EXPECT_NEAR(result_2, 1., 1e-7);
  EXPECT_TRUE(abserr_2 <= 1e-7);

  auto [result_3, abserr_3] = worker.qag_si(function1, -inf, inf, 1e-10, 1e-7);
  EXPECT_NEAR(result_3, 1.5, 1.5 * 1e-7);
  EXPECT_TRUE(abserr_3 <= 1.5 * 1e-7);

  auto [result_4, abserr_4] = worker.qag_si(function1, -1.0, 1.0, 1e-10, 1e-7);
  double result_4_analytic = 0.5 * (1. - std::exp(-2.)) + (1. - std::exp(-1.));
  EXPECT_NEAR(result_4, result_4_analytic, result_4_analytic * 1e-7);
  EXPECT_TRUE(abserr_4 <= result_4_analytic * 1e-7);

  gsl_set_error_handler(gsl_default_handler);
}

// *********

TEST(gsl_integration, complex_qag) { // NOLINT
  auto gsl_default_handler = gsl_set_error_handler_off();
  gsl_integration_wrapper worker(100);

  auto function = [](double t) -> dcomplex { return t * t + 1.0i * t; };
  auto [result, abserr] = worker.qag(function, 1, 2, 1e-10, 1e-7, GSL_INTEG_GAUSS31);
  EXPECT_NEAR(std::real(result), 7. / 3., 7. / 3. * 1e-7);
  EXPECT_NEAR(std::imag(result), 3. / 2., 3. / 2. * 1e-7);
  EXPECT_TRUE(std::real(abserr) < 7. / 3. * 1e-7);
  EXPECT_TRUE(std::imag(abserr) < 3. / 2. * 1e-7);

  gsl_set_error_handler(gsl_default_handler);
}

TEST(gsl_integration, complex_qag_si) { // NOLINT
  auto gsl_default_handler = gsl_set_error_handler_off();
  gsl_integration_wrapper worker(100);

  double inf = std::numeric_limits<double>::infinity();

  auto function = [](double t) -> dcomplex { return (int(t > 0.) + 1) * std::exp(-(1. + 1.0i) * t * t); };
  auto [result, abserr] = worker.qag_si(function, -inf, inf, 1e-10, 1e-7);

  dcomplex analytic_result = 3. / 2 * std::sqrt(M_PI / (1.0 + 1.0i));

  EXPECT_NEAR(std::real(result), std::real(analytic_result), std::abs(std::real(analytic_result)) * 1e-7);
  EXPECT_NEAR(std::imag(result), std::imag(analytic_result), std::abs(std::imag(analytic_result)) * 1e-7);
  EXPECT_TRUE(std::real(abserr) < std::abs(std::real(analytic_result)) * 1e-7);
  EXPECT_TRUE(std::imag(abserr) < std::abs(std::imag(analytic_result)) * 1e-7);

  gsl_set_error_handler(gsl_default_handler);
}

MAKE_MAIN // NOLINT
