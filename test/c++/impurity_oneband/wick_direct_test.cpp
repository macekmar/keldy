#include <keldy/common.hpp>
#include <keldy/impurity_oneband/wick_direct.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(integrand_direct, domain_checking) { // NOLINT
  model_param_t params;
  params.time_max = 20.0;
  params.Gamma = 0.5;
  params.beta = -1;
  params.nr_time_points_gf = 10001;
  params.bath_type = "flatband";
  params.ft_method = "analytic";

  g0_keldysh_contour_t g0_k{g0_model{g0_model_omega{params}, false}};
  auto external_A = gf_index_t{10.0, up, forward};
  auto external_B = gf_index_t{10.0, up, forward};

  integrand_g_direct integrand(g0_k, external_A, external_B, 0.); // no cutoff

  for (double t : {0., 0.5, 1., 5., 10.}) {
    auto [val, in_domain] = integrand(std::vector<double>{t, 3 * t, 2 * t});

    std::cout << "value: " << val << std::endl;
  }
}

MAKE_MAIN // NOLINT
