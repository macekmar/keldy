//#include <gtest/gtest.h>
#include <keldy/interfaces/gsl_qag_wrap.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy::details;

TEST(gsl_qag, real_integral) { // NOLINT

  gsl_qag_wrapper_t worker(100);

  auto function1 = [](double t) -> double { return t * t; };
  worker.integrate(function1, 1, 2, 1e-10, 1e-7, GSL_INTEG_GAUSS31);
  std::cout << "Result = " << worker.get_result() << " +- " << worker.get_abserr() << std::endl;
  EXPECT_NEAR(worker.get_result(), 2. + 1. / 3., 1e-6);

  auto function2 = [](double t) -> double { return 1. / (t * t); };
  worker.integrate(function2, 1, 2, 1e-10, 1e-7, GSL_INTEG_GAUSS31);
  std::cout << "Result = " << worker.get_result() << " +- " << worker.get_abserr() << std::endl;
  EXPECT_NEAR(worker.get_result(), 1. / 2., 1e-6);
}

TEST(gsl_qag, complex_integral) { // NOLINT

  gsl_qag_cpx_wrapper_t worker(100);

  auto function1 = [](double t) -> dcomplex { return t * t + 1_j * t; };
  worker.integrate(function1, 1, 2, 1e-10, 1e-7, GSL_INTEG_GAUSS31);
  std::cout << "Result = " << worker.get_result() << " +- " << worker.get_abserr_real() << ", "
            << worker.get_abserr_imag() << std::endl;
  EXPECT_COMPLEX_NEAR(worker.get_result(), 2. + 1. / 3. + 1_j * 3. / 2., 1e-6);

  auto function2 = [](double t) -> dcomplex { return 1. / (t * t) + 1_j * t; };
  worker.integrate(function2, 1, 2, 1e-10, 1e-7, GSL_INTEG_GAUSS31);
  std::cout << "Result = " << worker.get_result() << " +- " << worker.get_abserr_real() << ", "
            << worker.get_abserr_imag() << std::endl;
  EXPECT_COMPLEX_NEAR(worker.get_result(), 1. / 2. + 1_j * 3. / 2., 1e-6);
}

MAKE_MAIN; // NOLINT
