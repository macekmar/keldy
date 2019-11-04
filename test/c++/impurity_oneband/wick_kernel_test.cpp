#include "keldy/common.hpp"
#include <keldy/impurity_oneband/compute_obs.hpp>
#include <triqs/test_tools/gfs.hpp>
#include <algorithm>

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(integrand_kernel, Order_1) { // NOLINT
  model_param_t const params;
  g0_model const g0{params};
  g0_keldysh_contour_t const g0_k{g0};
  auto const external_B = gf_index_t{5.0, up, forward};
  integrand_g_kernel const integrand{g0_k, external_B};

  double const u = 2.0;

  auto expected_val =
     sparse_kernel_binner{
        {{{u, up, forward}, -g0_k({u, up, forward}, external_B) * g0_k({u, down, forward}, {u, down, backward})},
         {{u, up, backward}, g0_k({u, up, backward}, external_B) * g0_k({u, down, forward}, {u, down, backward})}}}
        .data;

  auto computed_val = integrand(std::vector<double>{u}).data;
  ASSERT_EQ(computed_val.size(), 2);

  auto const compare = [](std::pair<gf_index_t, dcomplex> const &a, std::pair<gf_index_t, dcomplex> const &b) -> bool {
    return a.first < b.first;
  };
  std::sort(expected_val.begin(), expected_val.end(), compare);
  std::sort(computed_val.begin(), computed_val.end(), compare);

  for (int i = 0; i < 2; ++i) {
    EXPECT_COMPLEX_NEAR(computed_val[i].second, expected_val[i].second, 1e-16);
  }
}

TEST(ComputeObs, Initialize2) { // NOLINT
  model_param_t params;
  compute_gf_kernel computer(params, 10.0, 4, "identity", 1000);

  computer.run(10);
  //

  // std::cout << computer.reduce_result() << std::endl;
  std::cout << computer.reduce_nr_points_run() << std::endl;
}

TEST(ComputeObs, Initialize3) { // NOLINT
  model_param_t params;
  compute_gf_kernel computer(params, 10.0, 4, "first_order", 1000);

  computer.run(10);

  auto result = computer.reduce_result();

  std::cout << result.get_bin_times() << std::endl;
  std::cout << result.get_bin_size() << std::endl;
  std::cout << result.get_nr_point_dropped() << std::endl;
  std::cout << result.get_nr_values()() << std::endl;
  std::cout << result.get_bin_times() << std::endl;

  //

  // std::cout << computer.reduce_result() << std::endl;
  std::cout << computer.reduce_nr_points_run() << std::endl;
}

MAKE_MAIN; // NOLINT
