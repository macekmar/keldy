#include "keldy/common.hpp"
#include <keldy/impurity_oneband/compute_obs.hpp>

#include <triqs/test_tools/gfs.hpp>
#include <algorithm>

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(integrand_kernel, Order_1) { // NOLINT
  model_param_t const params;
  g0_model const g0{g0_model_omega{params}, true};
  g0_keldysh_contour_t const g0_k{g0};
  auto const external_B = gf_index_t{5.0, up, forward};
  integrand_g_kernel const integrand{g0_k, external_B};

  double const u = 2.0;

  auto expected_val =
     sparse_kernel_binner{
        {{{u, up, forward}, -g0_k({u, up, forward}, external_B) * g0_k({u, down, forward}, {u, down, backward})},
         {{u, up, backward}, g0_k({u, up, backward}, external_B) * g0_k({u, down, forward}, {u, down, backward})}}}
        .data;

  auto computed_val = integrand(std::vector<double>{u}).first.data;
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

TEST(integrand_kernel, Order_2) { // NOLINT
  model_param_t const params;
  g0_model const g0{g0_model_omega{params}, true};
  g0_keldysh_contour_t const g0_k{g0};
  auto const external_B = gf_index_t{5.0, up, forward};
  integrand_g_kernel const integrand{g0_k, external_B};

  double const u = 2.0;
  double const v = 4.0;

  auto const det_down = [&g0_k, u, v](keldysh_idx_t a, keldysh_idx_t b) -> dcomplex {
    return g0_k({u, down, forward}, {u, down, backward}) * g0_k({v, down, forward}, {v, down, backward})
       - g0_k({v, down, b}, {u, down, a}) * g0_k({u, down, a}, {v, down, b});
  };

  auto const det_prod_u = [&g0_k, &det_down, u, v, external_B](keldysh_idx_t a, keldysh_idx_t b) -> dcomplex {
    int const sign = (a == b) ? +1 : -1;
    return sign
       * (g0_k({u, up, a}, {v, up, b}) * g0_k({v, up, b}, external_B) * det_down(a, b)
          + g0_k({v, up, forward}, {v, up, backward}) * g0_k({u, up, a}, external_B) * g0_k({v, down, b}, {u, down, a})
             * g0_k({u, down, a}, {v, down, b}));
  };

  auto const det_prod_v = [&g0_k, &det_down, u, v, external_B](keldysh_idx_t a, keldysh_idx_t b) -> dcomplex {
    int const sign = (a == b) ? +1 : -1;
    return -sign
       * (-g0_k({v, up, b}, {u, up, a}) * g0_k({u, up, a}, external_B) * det_down(a, b)
          - g0_k({u, up, forward}, {u, up, backward}) * g0_k({v, up, b}, external_B) * g0_k({v, down, b}, {u, down, a})
             * g0_k({u, down, a}, {v, down, b}));
  };

  auto expected_val =
     sparse_kernel_binner{{{{u, up, forward}, det_prod_u(forward, forward) + det_prod_u(forward, backward)},
                           {{u, up, backward}, det_prod_u(backward, forward) + det_prod_u(backward, backward)},
                           {{v, up, forward}, det_prod_v(forward, forward) + det_prod_v(backward, forward)},
                           {{v, up, backward}, det_prod_v(forward, backward) + det_prod_v(backward, backward)}}}
        .data;

  auto computed_val = integrand(std::vector<double>{u, v}).first.data;
  ASSERT_EQ(computed_val.size(), 4);

  auto const compare = [](std::pair<gf_index_t, dcomplex> const &a, std::pair<gf_index_t, dcomplex> const &b) -> bool {
    return a.first < b.first;
  };
  std::sort(expected_val.begin(), expected_val.end(), compare);
  std::sort(computed_val.begin(), computed_val.end(), compare);

  for (int i = 0; i < 4; ++i) {
    std::cout << computed_val[i].first << " ; " << expected_val[i].first << std::endl;
    std::cout << computed_val[i].second << " ; " << expected_val[i].second << std::endl;
    EXPECT_COMPLEX_NEAR(computed_val[i].second, expected_val[i].second, 1e-16);
  }
}

TEST(ComputeObs, Initialize2) { // NOLINT
  model_param_t params;
  double time = 10.0;
  compute_gf_kernel computer(params, time, 4);
  computer.warper.emplace_back(warpers::warper_plasma_uv_t{time});

  computer.run(10);
  //

  // std::cout << computer.reduce_result() << std::endl;
  std::cout << computer.reduce_nr_points_run() << std::endl;
}

TEST(ComputeObs, Initialize3) { // NOLINT

// todo: test something

  // model_param_t params;
  // compute_gf_kernel computer(params, 10.0, 4);
  // computer.warper.emplace_back(warpers::warper_plasma_uv_t{time});

  // computer.run(10);

  // auto result = computer.reduce_result();

  // std::cout << result.get_bin_times() << std::endl;
  // std::cout << result.get_bin_size() << std::endl;
  // std::cout << result.get_nr_point_dropped() << std::endl;
  // std::cout << result.get_nr_values()() << std::endl;
  // std::cout << result.get_bin_times() << std::endl;

  // //

  // // std::cout << computer.reduce_result() << std::endl;
  // std::cout << computer.reduce_nr_points_run() << std::endl;
}

MAKE_MAIN; // NOLINT
