#include <keldy/impurity_oneband/compute_obs.hpp>
#include <triqs/test_tools/gfs.hpp>
//#include <cmath>
#include <boost/math/special_functions/expint.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;

double const pi = 3.1415926535897932384626433832795028841971693993751058209749;

TEST(g0_model, Initialize) { // NOLINT
  model_param_t params;
  g0_model g0{params};
}

/*
 * Particle-hole symmetric case at T=0
 */
TEST(g0_model, Values) { // NOLINT

  model_param_t params;
  params.beta=100000.0;
  params.bias_V=0.0;
  params.eps_d=0.0;
  params.Gamma=1.0;
  params.time_max=1000.0;
  params.nr_time_points_gf=100000;
  params.alpha=0.0;
  params.bath_type="flatband";

  // TODO: copy params ?
  g0_model g0{params};

  auto g0_less_expected = [] (double t1, double t2) -> dcomplex {
    auto Gt = 1.0 * (t1 - t2); // = Gamma * (t1 - t2)
    auto real_part =
       (Gt == 0) ? 0.0 : (std::exp(Gt) * boost::math::expint(-Gt) - std::exp(-Gt) * boost::math::expint(Gt)) / (2 * pi);
    return real_part + 0.5_j * std::exp(-std::abs(Gt));
  };
  auto g0_grea_expected = [g0_less_expected] (double t1, double t2) -> dcomplex { return std::conj(g0_less_expected(t1, t2)); };

  double time;
  for (size_t i=0; i<=100; ++i) {
    time = -100.0 + i * 2.0;
    std::cout << "t = " << time << std::endl;
    EXPECT_COMPLEX_NEAR(g0_less_expected(time, 0.0), g0.g0_lesser[up](time), 1e-3);
    EXPECT_COMPLEX_NEAR(g0_grea_expected(time, 0.0), g0.g0_greater[up](time), 1e-3);
  }

}

TEST(g0_keldysh_adaptor, Initialize) { // NOLINT
  model_param_t params;
  g0_model g0{params};
  g0_keldysh_contour_t g0_k{g0};
}

MAKE_MAIN; // NOLINT
