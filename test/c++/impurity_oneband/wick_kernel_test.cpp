#include "keldy/common.hpp"
#include <keldy/impurity_oneband/compute_obs.hpp>

#include <triqs/test_tools/gfs.hpp>
#include <algorithm>

using namespace keldy;
using namespace keldy::impurity_oneband;

/// Sort the binner data in lexicographic order of the coordinates, so as to
//allow checking equality
template <int N, int M>
void sort(binner::sparse_binner_t<N, M> &sp_bin) {
  using data_t = std::pair<typename binner::sparse_binner_t<N, M>::coord_arr_t, dcomplex>;
  auto comp = [](data_t const &a, data_t const &b) -> bool {
    if (a.first.first == b.first.first) {
      return std::lexicographical_compare(a.first.second.cbegin(), a.first.second.cend(), b.first.second.cbegin(),
                                          b.first.second.cend());
    }
    return std::lexicographical_compare(a.first.first.cbegin(), a.first.first.cend(), b.first.first.cbegin(),
                                        b.first.first.cend());
  };
  std::stable_sort(sp_bin.data.begin(), sp_bin.data.end(), comp);
};

TEST(integrand_kernel, Order_1) { // NOLINT
  model_param_t const params;
  g0_model const g0{g0_model_omega{params}, true};
  g0_keldysh_contour_t const g0_k{g0};
  auto const external_B = gf_index_t{5.0, up, forward};
  integrand_g_kernel const integrand{g0_k, external_B};

  double const u = 2.0;

  auto expected_res = binner::sparse_binner_t<1, 1>();
  expected_res.accumulate(-g0_k({u, up, forward}, external_B) * g0_k({u, down, forward}, {u, down, backward}), u,
                          forward);
  expected_res.accumulate(g0_k({u, up, backward}, external_B) * g0_k({u, down, forward}, {u, down, backward}), u,
                          backward);
  sort(expected_res);

  binner::sparse_binner_t<1, 1> computed_res = integrand(std::vector<double>{u}).first;
  sort(computed_res);
  ASSERT_EQ(computed_res.data.size(), 2);

  for (int i = 0; i < 2; ++i) {
    EXPECT_COMPLEX_NEAR(computed_res.data[i].second, expected_res.data[i].second, 1e-16);
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

  auto expected_res = binner::sparse_binner_t<1, 1>();
  expected_res.accumulate(det_prod_u(forward, forward), u, forward);
  expected_res.accumulate(det_prod_u(forward, backward), u, forward);
  expected_res.accumulate(det_prod_u(backward, forward), u, backward);
  expected_res.accumulate(det_prod_u(backward, backward), u, backward);
  expected_res.accumulate(det_prod_v(forward, forward), v, forward);
  expected_res.accumulate(det_prod_v(backward, forward), v, forward);
  expected_res.accumulate(det_prod_v(forward, backward), v, backward);
  expected_res.accumulate(det_prod_v(backward, backward), v, backward);
  sort(expected_res);

  binner::sparse_binner_t<1, 1> computed_res = integrand(std::vector<double>{u, v}).first;
  sort(computed_res);

  ASSERT_EQ(computed_res.data.size(), 8);
  ASSERT_EQ(expected_res.data.size(), 8);

  for (int i = 0; i < 4; ++i) {
    std::cout << computed_res.data[2 * i].first << " ; " << expected_res.data[2 * i].first << std::endl;
    std::cout << computed_res.data[2 * i + 1].first << " ; " << expected_res.data[2 * i + 1].first << std::endl;
    std::cout << computed_res.data[2 * i].second << " ; " << expected_res.data[2 * i].second << std::endl;
    std::cout << computed_res.data[2 * i + 1].second << " ; " << expected_res.data[2 * i + 1].second << std::endl;
    EXPECT_COMPLEX_NEAR(computed_res.data[2 * i].second + computed_res.data[2 * i + 1].second,
                        expected_res.data[2 * i].second + expected_res.data[2 * i + 1].second, 1e-16);
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
