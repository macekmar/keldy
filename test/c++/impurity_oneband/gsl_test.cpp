#include <keldy/common.hpp>
#include <keldy/impurity_oneband/compute_obs.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(integrand_direct, Order_1) { // NOLINT
  model_param_t params;
  params.beta = 100000.0;
  params.bias_V = 0.0;
  params.eps_d = 1.0;
  params.Gamma = 1.0;
  params.time_max = 1000.0;
  params.nr_time_points_gf = 100000;
  params.alpha = 0.0;
  params.bath_type = "flatband";
  params.ft_method = "fft";

  compute_charge_Q_direct_gsl_vegas comp{params, 10.0, 3, "sobol"};

  comp.run(1000);

  std::cout << comp.get_result() << ", " << comp.get_error() << std::endl;

  // g0_model g0{params};
  // g0_keldysh_contour_t g0_k{g0};
  // auto external_A = gf_index_t{5.0, up, forward};
  // auto external_B = gf_index_t{5.0, up, forward};
  // integrand_g_direct integrand{g0_k, external_A, external_B};

  // dcomplex computed_val;
  // dcomplex expected_val;

  // // We ignore disconnected diagrams
  // expected_val = 0;

  // for (keldysh_idx_t a : {forward, backward}) {
  //   int sign = (a == 0) ? 1 : -1;
  //   expected_val += -sign * g0_k({1.0, down, forward}, {1.0, down, backward}) * g0_k(external_A, {1.0, up, a})
  //      * g0_k({1.0, up, a}, external_B);
  // }

  // computed_val = integrand(std::vector<double>{1.0});

  // EXPECT_COMPLEX_NEAR(expected_val, computed_val, 1e-17);
}
