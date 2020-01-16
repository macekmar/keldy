#include "keldy/common.hpp"
#include "keldy/impurity_oneband/model.hpp"
#include "keldy/contour_integral.hpp"
#include <gtest/gtest.h>
#include <itertools/itertools.hpp>
#include <keldy/impurity_oneband/compute_obs.hpp>
//#include <keldy/impurity_oneband/model.hpp>
#include <triqs/test_tools/gfs.hpp>
//#include <cmath>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/constants/constants.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(g0_model, Initialize_semicircle) { // NOLINT
  model_param_t params;
  params.bath_type = "semicircle";
  params.ft_method = "fft";

  g0_model g0{g0_model_omega{params}, true};
}

/*
 * Semi-circular epsilon_d=-1 at beta=10 and bias_V=2.0
 */
TEST(g0_model, Semicirc_fft) { // NOLINT

  model_param_t params;
  params.beta = 10.0;
  params.bias_V = 2.0;
  params.eps_d = -1.0;
  params.Gamma = 1.5;
  params.time_max = 100.0;
  params.nr_time_points_gf = 10000;
  params.alpha = 0.0;
  params.bath_type = "semicircle";
  params.ft_method = "fft";

  g0_model g0{g0_model_omega{params}, true};

  /// values from ctint_keldysh
  EXPECT_COMPLEX_NEAR(-0.23107349410223865 + 0.13213346415985422_j, g0.g0_lesser[up](1.0)(0, 0), 1e-3);
  EXPECT_COMPLEX_NEAR(0.17154101374667907 + 0.13375649387573293_j, g0.g0_lesser[up](1.0)(0, 1), 1e-3);
  EXPECT_COMPLEX_NEAR(-0.0066645142937036958 + 0.12215571972575601_j, g0.g0_lesser[up](1.0)(1, 0), 1e-3);
}

/*
 * Semi-circular epsilon_d=-1 at beta=10 and bias_V=2.0
 */
TEST(g0_model, Semicirc_contour) { // NOLINT

  model_param_t params;
  params.beta = 10.0;
  params.bias_V = 2.0;
  params.eps_d = -1.0;
  params.Gamma = 1.5;
  params.time_max = 100.0;
  params.nr_time_points_gf = 1000;
  params.alpha = 0.0;
  params.bath_type = "semicircle";
  params.ft_method = "contour";

  g0_model g0{g0_model_omega{params}, true};

  /// values from ctint_keldysh
  EXPECT_COMPLEX_NEAR(-0.23107349410223865 + 0.13213346415985422_j, g0.g0_lesser[up](1.0)(0, 0), 1e-3);
  //EXPECT_COMPLEX_NEAR(0.17154101374667907 + 0.13375649387573293_j, g0.g0_lesser[up](1.0)(0, 1), 1e-3);
  EXPECT_COMPLEX_NEAR(-0.0066645142937036958 + 0.12215571972575601_j, g0.g0_lesser[up](1.0)(1, 0), 1e-3);
}

/// Difference between consecutive elements of a 1D array
template <typename T>
auto diff(T const &arr) {
  auto result = array<typename T::value_type, 1>(first_dim(arr) - 1);
  foreach (result, [&arr, &result](size_t i) { result(i) = arr(i + 1) - arr(i); })
    ;
  return result;
}

/*
 * Test presence of outliers
 */
TEST(g0_model, Semicirc_contour_outliers) { // NOLINT

  model_param_t params;
  params.beta = -1.;
  params.bias_V = -3.0;
  params.eps_d = -0.5;
  params.Gamma = 1.;
  params.time_max = 50.0;
  params.nr_time_points_gf = 5001;
  params.alpha = 0.0;
  params.bath_type = "semicircle";
  params.ft_method = "contour";

  g0_model g0{g0_model_omega{params}, true};

  // With these parameters, the second derivative of the on-site g0 should be lower than 1 in amplitude
  auto data = g0.g0_lesser[up].data()(range(), 0, 0);
  double const dt = g0.g0_lesser[up].mesh().delta();
  auto deriv = diff(data) / dt;
  auto deriv2 = diff(deriv) / dt;

  EXPECT_LT(max_element(abs(real(deriv2))), 1.);
  EXPECT_LT(max_element(abs(imag(deriv2))), 1.);
}

MAKE_MAIN; // NOLINT
