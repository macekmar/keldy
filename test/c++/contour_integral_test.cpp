#include "keldy/common.hpp"
#include "keldy/contour_integral.hpp"
#include <gtest/gtest.h>
#include <itertools/itertools.hpp>
#include <triqs/test_tools/gfs.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/constants/constants.hpp>

using namespace keldy;
using namespace std::complex_literals;

TEST(g0_model, contour_integration) { // NOLINT
  using namespace boost::math::double_constants;
  auto gsl_error_handler_old = gsl_set_error_handler_off();
  contour_integration_t worker_1(-3., 2.);

  auto f = [](dcomplex x) -> dcomplex { return std::exp(-(3. + 1.0i) * x * x); };
  dcomplex const ref = 1. / std::sqrt((3. + 1.0i) / pi);

  worker_1.integrate(f, -1., +1.);
  EXPECT_COMPLEX_NEAR(worker_1.get_result(), ref, 1e-8);

  worker_1.integrate(f, -2. - 1.i, +2. - 1.i);
  EXPECT_COMPLEX_NEAR(worker_1.get_result(), ref, 1e-8);

  contour_integration_t worker_2(-3., 2.);

  worker_2.integrate(f, -1., +1.);
  EXPECT_COMPLEX_NEAR(worker_2.get_result(), ref, 1e-8);

  worker_2.integrate(f, -2. - 1.i, +2. - 1.i);
  EXPECT_COMPLEX_NEAR(worker_2.get_result(), ref, 1e-8);

  gsl_set_error_handler(gsl_error_handler_old);
}

MAKE_MAIN // NOLINT
