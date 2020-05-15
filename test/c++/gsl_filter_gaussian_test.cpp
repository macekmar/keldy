#include <gtest/gtest.h>
#include <keldy/interfaces/gsl_filter_gaussian_wrap.hpp>
#include <triqs/test_tools/arrays.hpp>

using namespace triqs::arrays;
using namespace keldy::details;

TEST(GSLFilterGaussian, Value) { //NOLINT
  array<double, 1> data(100);
  for (int i = 0; i < data.size(); ++i) {
    data(i) = -0.04 * i * i + 3.0 * i + 6.0;
  }
  double const sigma = 1.5;

  auto filtered_data = gsl_filter_gaussian_wrapper(data, 10, sigma, GSL_FILTER_END_PADZERO);
  EXPECT_EQ(filtered_data.size(), 100);

  auto kernel = gsl_filter_gaussian_kernel_wrapper(10, sigma, false);
  double norm = sum(kernel);

  EXPECT_DOUBLE_EQ(filtered_data(65), sum(data(range(60, 71)) * kernel) / norm);
  EXPECT_DOUBLE_EQ(filtered_data(98), sum(data(range(93, 100)) * kernel(range(0, 7))) / norm);
}

TEST(GSLFilterGaussian, Kernel) { //NOLINT
  auto kernel = gsl_filter_gaussian_kernel_wrapper(40, 4.5, false);
  EXPECT_EQ(kernel.size(), 41);
  EXPECT_DOUBLE_EQ(kernel(20), 1.0); // not normalized
}

MAKE_MAIN; // NOLINT
