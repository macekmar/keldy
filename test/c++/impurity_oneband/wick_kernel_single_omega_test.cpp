#include <keldy/common.hpp>
#include <keldy/impurity_oneband/wick_direct.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(kernel_single_omega, integrand) { // NOLINT
  model_param_t params;
  params.beta = 1.0;
  params.bias_V = 0.0;
  params.eps_d = 0.0;
  params.Gamma = 1.0;
  params.time_max = 10.0; // (time_limit_min = - time_limit_max)
  params.nr_time_points_gf = 1000;
  params.alpha = 0.0;
  params.bath_type = "flatband";
  params.ft_method = "fft";

  g0_model g0{g0_model_omega{params}, false};
  g0_keldysh_contour_t g0_k{g0};
  auto external_B = gf_index_t{5.0, up, forward};
  double omega = 1.;

  integrand_g_kernel_single_omega integrand(g0_k, external_B, omega);

  std::vector<double> time_vec = {2., 0.};
  std::cout << integrand(time_vec) << std::endl;

  time_vec = {4.9, 4.8};
  std::cout << integrand(time_vec) << std::endl;
}

MAKE_MAIN // NOLINT
