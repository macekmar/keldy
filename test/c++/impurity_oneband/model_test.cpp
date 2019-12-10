#include "keldy/common.hpp"
#include "keldy/impurity_oneband/model.hpp"
#include <gtest/gtest.h>
#include <itertools/itertools.hpp>
#include <keldy/impurity_oneband/compute_obs.hpp>
#include <triqs/test_tools/gfs.hpp>
//#include <cmath>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/constants/constants.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(g0_model, Initialize_flatband) { // NOLINT
  model_param_t params;
  params.bath_type = "flatband";
  g0_model g0{params, true};
}

TEST(g0_model, Initialize_flatband_analytic) { // NOLINT
  model_param_t params;
  params.beta = -1.;
  params.bias_V = 0.;
  params.eps_d = 0.;
  params.alpha = 0.;
  params.bath_type = "flatband_analytic";
  g0_model g0{params, false};
}

TEST(g0_model, Initialize_semicircle) { // NOLINT
  model_param_t params;
  params.bath_type = "semicircle";
  g0_model g0{params, true};
}

/*
 * Flat band particle-hole symmetric case at T=0
 */
TEST(g0_model, Flatband_sym) { // NOLINT

  model_param_t params;
  params.beta = 100000.0;
  params.bias_V = 0.0;
  params.eps_d = 0.0;
  params.Gamma = 1.5;
  params.time_max = 1000.0;
  params.nr_time_points_gf = 100000;
  params.alpha = 0.0;
  params.bath_type = "flatband";

  // TODO: copy params ?
  g0_model g0{params, false};

  auto g0_less_expected = [Gamma = params.Gamma](double t1, double t2) -> dcomplex {
    using namespace boost::math::double_constants;
    auto Gt = Gamma * (t1 - t2); // = Gamma * (t1 - t2)
    auto real_part =
       (Gt == 0) ? 0.0 : (std::exp(Gt) * boost::math::expint(-Gt) - std::exp(-Gt) * boost::math::expint(Gt)) / (2 * pi);
    return real_part + 0.5_j * std::exp(-std::abs(Gt));
  };
  auto g0_grea_expected = [g0_less_expected](double t1, double t2) -> dcomplex {
    return std::conj(g0_less_expected(t1, t2));
  };

  double time;
  for (size_t i = 0; i <= 100; ++i) {
    time = -100.0 + i * 2.0;
    std::cout << "t = " << time << std::endl;
    EXPECT_COMPLEX_NEAR(g0_less_expected(time, 0.0), g0.g0_lesser[up](time)(0, 0), 1e-3);
    EXPECT_COMPLEX_NEAR(g0_grea_expected(time, 0.0), g0.g0_greater[up](time)(0, 0), 1e-3);
  }
}

/*
 * Flat band epsilon_d=1 at T=0
 */
TEST(g0_model, Flatband_asym_1) { // NOLINT

  model_param_t params;
  params.beta = 100000.0;
  params.bias_V = 0.0;
  params.eps_d = 1.0;
  params.Gamma = 1.5;
  params.time_max = 100.0;
  params.nr_time_points_gf = 10000;
  params.alpha = 0.0;
  params.bath_type = "flatband";

  g0_model g0{params, false};

  /// value from Marjan analytic formula
  EXPECT_COMPLEX_NEAR(-0.09907315 + 0.06406786_j, g0.g0_lesser[up](1.0)(0, 0), 1e-4);
}

/*
 * Flat band epsilon_d=-1 at T=0
 */
TEST(g0_model, Flatband_asym_2) { // NOLINT

  model_param_t params;
  params.beta = 100000.0;
  params.bias_V = 0.0;
  params.eps_d = -1.0;
  params.Gamma = 1.5;
  params.time_max = 100.0;
  params.nr_time_points_gf = 10000;
  params.alpha = 0.0;
  params.bath_type = "flatband";

  g0_model g0{params, false};

  /// value from Marjan analytic formula
  EXPECT_COMPLEX_NEAR(-0.2868307 + 0.05648988_j, g0.g0_lesser[up](1.0)(0, 0), 1e-4);
}

/*
 * Flat band epsilon_d=-1 at beta=10 and bias_V=2.0
 */
TEST(g0_model, Flatband_3) { // NOLINT

  model_param_t params;
  params.beta = 10.0;
  params.bias_V = 2.0;
  params.eps_d = -1.0;
  params.Gamma = 1.5;
  params.time_max = 100.0;
  params.nr_time_points_gf = 10000;
  params.alpha = 0.0;
  params.bath_type = "flatband";

  g0_model g0{params, true};

  /// values from ctint_keldysh (dot-lead terms have an extra factor 1/2, mistake in the original code)
  EXPECT_COMPLEX_NEAR(-0.218086 + 0.0265308_j, g0.g0_lesser[up](1.0)(0, 0), 1e-3);
  EXPECT_COMPLEX_NEAR((0.258951 + 0.423124_j) / 2, g0.g0_lesser[up](1.0)(0, 1), 1e-3);
  EXPECT_COMPLEX_NEAR((-0.107289 + 0.350721_j) / 2, g0.g0_lesser[up](1.0)(1, 0), 1e-3);
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

  g0_model g0{params, true};

  /// values from ctint_keldysh
  EXPECT_COMPLEX_NEAR(-0.23107349410223865 + 0.13213346415985422_j, g0.g0_lesser[up](1.0)(0, 0), 1e-3);
  EXPECT_COMPLEX_NEAR(0.17154101374667907 + 0.13375649387573293_j, g0.g0_lesser[up](1.0)(0, 1), 1e-4);
  EXPECT_COMPLEX_NEAR(-0.0066645142937036958 + 0.12215571972575601_j, g0.g0_lesser[up](1.0)(1, 0), 1e-4);
}

TEST(g0_keldysh_adaptor, Initialize) { // NOLINT
  model_param_t params;
  g0_model g0{params, true};
  g0_keldysh_contour_t g0_k{g0};
}

MAKE_MAIN; // NOLINT
