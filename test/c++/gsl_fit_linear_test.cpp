#include <triqs/arrays.hpp>
#include <keldy/interfaces/gsl_fit_linear_wrap.hpp>
#include <gtest/gtest.h>
#include <triqs/test_tools/arrays.hpp>

using namespace triqs::arrays;
using namespace keldy::details;

TEST(GSLFitLinear, SimpleUsage) { // NOLINT
  array<double, 1> x_arr(10);
  array<double, 1> y_arr(10);
  array<double, 1> w_arr(10);

  for (auto [i, x] : itertools::enumerate(x_arr)) {
    x_arr(i) = i / 9.0 - 6;
    y_arr(i) = 2 * x + 3;
    w_arr(i) = 1 / (5 * 5 + (x - 3) * (x - 3));
  }

  auto res = gsl_fit_wlinear_wrapper(x_arr, y_arr, w_arr);

  std::cout << res.c0 << std::endl;
  std::cout << res.c1 << std::endl;
  std::cout << res.cov << std::endl;
  std::cout << res.chisq << std::endl;

  EXPECT_DOUBLE_EQ(res.c0, 3);
  EXPECT_DOUBLE_EQ(res.c1, 2);
  EXPECT_DOUBLE_EQ(res.chisq, 0);
}

MAKE_MAIN; // NOLINT
