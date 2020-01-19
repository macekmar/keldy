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
TEST(g0_model, Semicirc) { // NOLINT

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

MAKE_MAIN; // NOLINT
