#include "keldy/warpers/plasma_uv.hpp"
#include "keldy/warpers/product_1d_simple.hpp"
#include <gtest/gtest.h>
#include <keldy/common.hpp>
#include <keldy/impurity_oneband/compute_obs.hpp>
#include <keldy/impurity_oneband/wick_direct_time.hpp>
#include <keldy/impurity_oneband/wick_direct.hpp>
#include <triqs/test_tools/arrays.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;

// Complex Double equal
#define EXPECT_DCOMPLEX_COORD_EQ(a, b)                                                                                 \
  EXPECT_DOUBLE_EQ((a).real(), (b).real());                                                                            \
  EXPECT_DOUBLE_EQ((a).imag(), (b).imag())

TEST(integrand_direct_time, time) { // NOLINT
  model_param_t params;
  params.time_max = 10.0;
  params.nr_time_points_gf = 100;
  params.bath_type = "flatband";
  params.ft_method = "fft";

  g0_keldysh_contour_t g0_k{g0_model{g0_model_omega{params}, false}};
  auto external_A = gf_index_t{5.0, up, forward};
  auto external_B = gf_index_t{5.0, up, forward};

  integrand_g_direct_time integrand(g0_k, external_A, external_B, 0.); // no cutoff

  std::vector<double> time_vec = {2.5, 1.2};
  binner::sparse_binner_t<1> res = integrand(time_vec).first;
  EXPECT_EQ(res.data.size(), 1);
  EXPECT_EQ(res.data[0].first.first[0], 1.2);
}

TEST(integrand_direct_time, consistency) { // NOLINT
  model_param_t params;
  params.time_max = 10.0;
  params.nr_time_points_gf = 100;
  params.bath_type = "flatband";
  params.ft_method = "fft";

  g0_keldysh_contour_t g0_k{g0_model{g0_model_omega{params}, false}};
  auto external_A = gf_index_t{5.0, up, forward};
  auto external_B = gf_index_t{5.0, up, forward};

  integrand_g_direct_time integrand(g0_k, external_A, external_B, 0.);   // no cutoff
  integrand_g_direct integrand_direct(g0_k, external_A, external_B, 0.); // no cutoff

  std::vector<double> time_vec = {2.5, 1.2};
  binner::sparse_binner_t<1> res = integrand(time_vec).first;
  EXPECT_EQ(res.data.size(), 1);
  auto integrand_direct_out = integrand_direct(time_vec).first;

  EXPECT_COMPLEX_NEAR(res.data[0].second, integrand_direct_out, 1e-16);
}

TEST(integration_direct_time, consistency) { // NOLINT
  model_param_t params;
  double time = 10.0;
  compute_charge_Q_direct computer_direct(params, time, 4, 0.0);
  compute_charge_Q_direct_time computer(params, time, 4, 10, 0.0);

  computer_direct.warper.emplace_back(warpers::warper_plasma_uv_t{time});
  computer_direct.warper.emplace_back(warpers::make_product_1d_simple_exponential(time, 1.0));

  computer.warper.emplace_back(warpers::warper_plasma_uv_t{time});
  computer.warper.emplace_back(warpers::make_product_1d_simple_exponential(time, 1.0));

  computer_direct.run(10);
  computer.run(10);

  dcomplex result_direct = computer_direct.reduce_result();
  binner::binner_t<1> result = computer.reduce_result();

  std::cout << result_direct << std::endl;
  std::cout << result.get_data() << std::endl;

  auto integrand_summed = sum(result.get_data());

  EXPECT_COMPLEX_NEAR(integrand_summed, result_direct, 1e-16);
}

MAKE_MAIN; // NOLINT
