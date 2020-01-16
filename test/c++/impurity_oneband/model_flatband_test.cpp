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

TEST(g0_model, Initialize_flatband) { // NOLINT
  model_param_t params;
  params.bath_type = "flatband";
  params.ft_method = "fft";
  g0_model g0{g0_model_omega{params}, false};
}

TEST(g0_model, Initialize_flatband_contour) { // NOLINT
  model_param_t params;
  params.beta = 1.0;
  params.bias_V = 0.3;
  params.eps_d = -1.0;
  params.Gamma = 1.5;
  params.time_max = 100.0;
  params.nr_time_points_gf = 100;
  params.alpha = 0.0;
  params.bath_type = "flatband";
  params.ft_method = "contour";

  g0_model g0{g0_model_omega{params}, false};

  // large beta
  params.beta = 10000.0;
  g0_model g0_cold{g0_model_omega{params}, false};
}

TEST(g0_model, Initialize_flatband_analytic) { // NOLINT
  model_param_t params;
  params.beta = -1.;
  params.bias_V = 0.;
  params.eps_d = 0.;
  params.alpha = 0.;
  params.bath_type = "flatband";
  params.ft_method = "analytic";

  g0_model g0{g0_model_omega{params}, false};
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
  params.nr_time_points_gf = 100001;
  params.alpha = 0.0;
  params.bath_type = "flatband";
  params.ft_method = "fft";

  // TODO: copy params ?
  g0_model g0{g0_model_omega{params}, false};

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
 * Flat band particle-hole symmetric case at T=0
 * contour integration
 */
TEST(g0_model, Flatband_contour_sym) { // NOLINT

  model_param_t params;
  params.beta = -1.;
  params.bias_V = 0.0;
  params.eps_d = 0.0;
  params.Gamma = 1.5;
  params.time_max = 20.0;
  params.nr_time_points_gf = 1001;
  params.alpha = 0.0;
  params.bath_type = "flatband";
  params.ft_method = "contour";

  g0_model g0_zeroT{g0_model_omega{params}, false};

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
    time = -10.0 + i * 0.2;
    //std::cout << "t = " << time << std::endl;
    EXPECT_COMPLEX_NEAR(g0_less_expected(time, 0.0), g0_zeroT.g0_lesser[up](time)(0, 0), 1e-8);
    EXPECT_COMPLEX_NEAR(g0_grea_expected(time, 0.0), g0_zeroT.g0_greater[up](time)(0, 0), 1e-8);
  }

  std::cout << "error = " << g0_zeroT.lesser_ft_error << std::endl;
  std::cout << "error = " << g0_zeroT.greater_ft_error << std::endl;

  params.beta = 100000.;
  g0_model g0_largebeta{g0_model_omega{params}, false};

  for (size_t i = 0; i <= 100; ++i) {
    time = -10.0 + i * 0.2;
    //std::cout << "t = " << time << std::endl;
    EXPECT_COMPLEX_NEAR(g0_less_expected(time, 0.0), g0_largebeta.g0_lesser[up](time)(0, 0), 1e-8);
    EXPECT_COMPLEX_NEAR(g0_grea_expected(time, 0.0), g0_largebeta.g0_greater[up](time)(0, 0), 1e-8);
  }

  std::cout << "error = " << g0_largebeta.lesser_ft_error << std::endl;
  std::cout << "error = " << g0_largebeta.greater_ft_error << std::endl;
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
  params.nr_time_points_gf = 10001;
  params.alpha = 0.0;
  params.bath_type = "flatband";
  params.ft_method = "fft";

  g0_model g0_fft{g0_model_omega{params}, false};

  /// value from Marjan analytic formula
  EXPECT_COMPLEX_NEAR(-0.09907315 + 0.06406786_j, g0_fft.g0_lesser[up](1.0)(0, 0), 1e-4);

  params.beta = -1;
  params.time_max = 10.0;
  params.nr_time_points_gf = 10001;
  params.bath_type = "flatband";
  params.ft_method = "contour";

  g0_model g0_ctr{g0_model_omega{params}, false};

  /// value from Marjan analytic formula
  EXPECT_COMPLEX_NEAR(-0.09907315 + 0.06406786_j, g0_ctr.g0_lesser[up](1.0)(0, 0), 1e-7);
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
  params.nr_time_points_gf = 10001;
  params.alpha = 0.0;
  params.bath_type = "flatband";
  params.ft_method = "fft";

  g0_model g0_fft{g0_model_omega{params}, false};

  /// value from Marjan analytic formula
  EXPECT_COMPLEX_NEAR(-0.2868307 + 0.05648988_j, g0_fft.g0_lesser[up](1.0)(0, 0), 1e-4);

  params.beta = -1;
  params.time_max = 10.0;
  params.nr_time_points_gf = 10001;
  params.bath_type = "flatband";
  params.ft_method = "contour";

  g0_model g0_ctr{g0_model_omega{params}, false};

  /// value from Marjan analytic formula
  EXPECT_COMPLEX_NEAR(-0.2868307 + 0.05648988_j, g0_ctr.g0_lesser[up](1.0)(0, 0), 1e-7);
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
  params.nr_time_points_gf = 10001;
  params.alpha = 0.0;
  params.bath_type = "flatband";
  params.ft_method = "fft";

  g0_model g0_fft{g0_model_omega{params}, false};

  std::cout << "t=0: " << g0_fft.g0_lesser[up](0.0) << std::endl;
  std::cout << "t=0: " << g0_fft.g0_lesser[up](0.0)(1, 0) << std::endl;
  std::cout << "t=0.1: " << g0_fft.g0_lesser[up](0.1) << std::endl;

  /// values from ctint_keldysh (dot-lead terms have an extra factor 1/2, mistake in the original code)
  EXPECT_COMPLEX_NEAR(-0.218086 + 0.0265308_j, g0_fft.g0_lesser[up](1.0)(0, 0), 1e-3);
  //EXPECT_COMPLEX_NEAR((0.258951 + 0.423124_j) / 2, g0_fft.g0_lesser[up](1.0)(0, 1), 1e-3);
  //EXPECT_COMPLEX_NEAR((-0.107289 + 0.350721_j) / 2, g0_fft.g0_lesser[up](1.0)(1, 0), 1e-3);

  /// check t <-> -t symmetry
  EXPECT_COMPLEX_NEAR(g0_fft.g0_lesser[up](1.0)(0, 0), -std::conj(g0_fft.g0_lesser[up](-1.0)(0, 0)), 1e-8);
  //EXPECT_COMPLEX_NEAR(g0_fft.g0_lesser[up](1.0)(0, 1), -std::conj(g0_fft.g0_lesser[up](-1.0)(1, 0)), 1e-8);
  EXPECT_COMPLEX_NEAR(g0_fft.g0_greater[up](1.0)(0, 0), -std::conj(g0_fft.g0_greater[up](-1.0)(0, 0)), 1e-8);
  //EXPECT_COMPLEX_NEAR(g0_fft.g0_greater[up](1.0)(0, 1), -std::conj(g0_fft.g0_greater[up](-1.0)(1, 0)), 1e-8);

  params.bath_type = "flatband";
  params.ft_method = "contour";

  g0_model g0_ctr{g0_model_omega{params}, false};

  std::cout << "t=0: " << g0_ctr.g0_lesser[up](0.0) << std::endl;
  std::cout << "t=0: " << g0_ctr.g0_lesser[up](0.0)(1, 0) << std::endl;
  std::cout << "t=0.1: " << g0_ctr.g0_lesser[up](0.1) << std::endl;

  /// values from ctint_keldysh (dot-lead terms have an extra factor 1/2, mistake in the original code)
  EXPECT_COMPLEX_NEAR(-0.218086 + 0.0265308_j, g0_ctr.g0_lesser[up](1.0)(0, 0), 1e-3);
  //EXPECT_COMPLEX_NEAR((0.258951 + 0.423124_j) / 2, g0_ctr.g0_lesser[up](1.0)(0, 1), 1e-3);
  //EXPECT_COMPLEX_NEAR((-0.107289 + 0.350721_j) / 2, g0_ctr.g0_lesser[up](1.0)(1, 0), 1e-3);

  /// check t <-> -t symmetry
  EXPECT_COMPLEX_NEAR(g0_ctr.g0_lesser[up](1.0)(0, 0), -std::conj(g0_ctr.g0_lesser[up](-1.0)(0, 0)), 1e-8);
  //EXPECT_COMPLEX_NEAR(g0_ctr.g0_lesser[up](1.0)(0, 1), -std::conj(g0_ctr.g0_lesser[up](-1.0)(1, 0)), 1e-8);
  EXPECT_COMPLEX_NEAR(g0_ctr.g0_greater[up](1.0)(0, 0), -std::conj(g0_ctr.g0_greater[up](-1.0)(0, 0)), 1e-8);
  //EXPECT_COMPLEX_NEAR(g0_ctr.g0_greater[up](1.0)(0, 1), -std::conj(g0_ctr.g0_greater[up](-1.0)(1, 0)), 1e-8);
}

MAKE_MAIN; // NOLINT
