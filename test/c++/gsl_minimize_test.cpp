#include <gtest/gtest.h>
#include <keldy/interfaces/gsl_minimize.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy::details;

/// Make sure we reset GSL error handler for each test.
class GSLMinimize : public ::testing::Test {
 private:
  gsl_error_handler_t *previous_handler = NULL;

 protected:
  void SetUp() override {
    previous_handler = gsl_set_error_handler(NULL); // set default error handler
  }
  void TearDown() override { gsl_set_error_handler(previous_handler); }
};

TEST_F(GSLMinimize, minimize) { // NOLINT
  auto f = [](double x) { return x * x + 3 * x * x * x * x; };

  auto res = gsl_minimize(f, -2., 1., 2., 1e-10, 0., 100, true);

  std::cout << res;

  EXPECT_NEAR(res.x, 0., 1e-10);
  EXPECT_LE(res.nr_iter, 100);
  EXPECT_EQ(res.status, GSL_SUCCESS);
  EXPECT_DOUBLE_EQ(res.f, f(res.x));
  EXPECT_DOUBLE_EQ(res.f_lower, f(res.x_lower));
  EXPECT_DOUBLE_EQ(res.f_upper, f(res.x_upper));
  EXPECT_LT(res.x_lower, res.x);
  EXPECT_LT(res.x, res.x_upper);
  EXPECT_LT(res.f, res.f_lower);
  EXPECT_LT(res.f, res.f_upper);
}

TEST_F(GSLMinimize, reach_max_iter) { // NOLINT
  auto f = [](double x) { return x * x + 3 * x * x * x * x; };

  auto res = gsl_minimize(f, -2., 1.999, 2., 1e-10, 0., 2, true);

  std::cout << res;

  EXPECT_EQ(res.nr_iter, 2);
  EXPECT_EQ(res.status, GSL_CONTINUE);
  EXPECT_DOUBLE_EQ(res.f, f(res.x));
  EXPECT_DOUBLE_EQ(res.f_lower, f(res.x_lower));
  EXPECT_DOUBLE_EQ(res.f_upper, f(res.x_upper));
  EXPECT_LT(res.x_lower, res.x);
  EXPECT_LT(res.x, res.x_upper);
  EXPECT_LT(res.f, res.f_lower);
  EXPECT_LT(res.f, res.f_upper);
}

/// make sure we throw python-friendly triqs runtime error rather than gsl errors
TEST_F(GSLMinimize, no_minimum) { // NOLINT
  auto f = [](double x) { return x; };

  EXPECT_THROW(gsl_minimize(f, -2., 1., 2., 1e-10, 0., 100), triqs::runtime_error);
}

/// make sure we throw python-friendly triqs runtime error rather than gsl errors
TEST_F(GSLMinimize, invalid_function) { // NOLINT
  auto f = [](double x) {
    if (x <= 0 || x >= 1) {
      return x * x + 3 * x * x * x * x;
    }
    return std::nan("1");
  };

  EXPECT_THROW(gsl_minimize(f, -2., 1.5, 2., 1e-10, 0., 100), triqs::runtime_error);
}

/// make sure we throw python-friendly triqs runtime error rather than gsl errors
TEST_F(GSLMinimize, invalid_function_start) { // NOLINT
  auto f = [](double x) {
    if (x <= 0) {
      return x * x + 3 * x * x * x * x;
    }
    return std::nan("1");
  };

  EXPECT_THROW(gsl_minimize(f, -2., 1.5, 2., 1e-10, 0., 100), triqs::runtime_error);
}

MAKE_MAIN // NOLINT
