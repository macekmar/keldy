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

TEST(LocalLinearRegression, NormalKernel) { //NOLINT
  /// kernel.size() < y.size()
  array<double, 1> x(100);
  array<double, 1> y(100);
  array<double, 1> y_out(100);
  for (int i = 0; i < x.size(); ++i) {
    x(i) = i;
    y(i) = 3.0 * i + 6.0;
    y_out(i) = 0;
  }
  double const sigma = 15;
  array<double, 1> kernel(51);
  for (int j = 0; j < kernel.size(); ++j) {
    kernel(j) = 1. / ((j - 25) * (j - 25) + sigma * sigma);
  }

  local_linear_reg(x, y, y_out, kernel);

  EXPECT_ARRAY_NEAR(y, y_out);
}

TEST(LocalLinearRegression, WideKernel) { //NOLINT
  /// y.size() < kernel.size() < 2 * y.size()
  array<double, 1> x(100);
  array<double, 1> y(100);
  array<double, 1> y_out(100);
  for (int i = 0; i < x.size(); ++i) {
    x(i) = i;
    y(i) = 3.0 * i + 6.0;
    y_out(i) = 0;
  }
  double const sigma = 15;
  array<double, 1> kernel(151);
  for (int j = 0; j < kernel.size(); ++j) {
    kernel(j) = 1. / ((j - 75) * (j - 75) + sigma * sigma);
  }

  local_linear_reg(x, y, y_out, kernel);

  EXPECT_ARRAY_NEAR(y, y_out);
}

TEST(LocalLinearRegression, VeryWideKernel) { //NOLINT
  /// kernel.size() > 2 * y.size()
  array<double, 1> x(100);
  array<double, 1> y(100);
  array<double, 1> y_out(100);
  for (int i = 0; i < x.size(); ++i) {
    x(i) = i;
    y(i) = 3.0 * i + 6.0;
    y_out(i) = 0;
  }
  double const sigma = 15;
  array<double, 1> kernel(301);
  for (int j = 0; j < kernel.size(); ++j) {
    kernel(j) = 1. / ((j - 150) * (j - 150) + sigma * sigma);
  }

  local_linear_reg(x, y, y_out, kernel);

  EXPECT_ARRAY_NEAR(y, y_out);
}

TEST(LocalLinearRegression, DeadPoints) { //NOLINT
  array<double, 1> x(100);
  array<double, 1> y(100);
  array<double, 1> y_out(100);
  for (int i = 0; i < x.size(); ++i) {
    x(i) = i;
    y(i) = 3.0 * i + 6.0;
    y_out(i) = 0;
  }
  std::vector<int> dead_points = {0, 15, 50, 51, 52, 60};
  double const sigma = 15;
  array<double, 1> kernel(51);
  for (int j = 0; j < kernel.size(); ++j) {
    kernel(j) = 1. / ((j - 25) * (j - 25) + sigma * sigma);
  }

  local_linear_reg(x, y, y_out, kernel, dead_points);

  EXPECT_ARRAY_NEAR(y, y_out);
}

TEST(LocalLinearRegression, ConsistencyKernelSizes) { //NOLINT
  array<double, 1> x(100);
  array<double, 1> y(100);
  array<double, 1> y_out(100);
  array<double, 1> y_out2(100);
  for (int i = 0; i < x.size(); ++i) {
    x(i) = i;
    y(i) = 0.1 * i * i - 3.0 * i + 6.0;
    y_out(i) = 0;
    y_out2(i) = 0;
  }
  double const sigma = 50;
  array<double, 1> kernel2(251);
  for (int j = 0; j < kernel2.size(); ++j) {
    if (j < 100 or 151 <= j) {
      kernel2(j) = 0;
    } else {
      kernel2(j) = 1. / ((j - 125) * (j - 125) + sigma * sigma);
    }
  }
  array<double, 1> kernel(51);
  for (int k = 0; k < kernel.size(); ++k) {
    kernel(k) = 1. / ((k - 25) * (k - 25) + sigma * sigma);
  }

  EXPECT_EQ(sum(kernel), sum(kernel2));
  local_linear_reg(x, y, y_out, kernel);
  local_linear_reg(x, y, y_out2, kernel2);

  EXPECT_ARRAY_NEAR(y_out, y_out2);
}

///TODO: add test against values
///TODO: test dead points

MAKE_MAIN // NOLINT
