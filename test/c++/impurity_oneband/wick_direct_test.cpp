#include <keldy/common.hpp>
#include <keldy/impurity_oneband/wick_direct.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(integrand_direct, Order_1) { // NOLINT
  model_param_t params;
  g0_model g0{params};
  g0_keldysh_contour_t g0_k{g0};
  auto external_A = gf_index_t{5.0, up, forward};
  auto external_B = gf_index_t{5.0, up, forward};
  integrand_g_direct integrand{g0_k, external_A, external_B};

  dcomplex computed_val;
  dcomplex expected_val;

  // We ignore disconnected diagrams
  expected_val = 0;

  for (keldysh_idx_t a : {forward, backward}) {
    int sign = (a == 0) ? 1 : -1;
    expected_val += -sign * g0_k({1.0, down, forward}, {1.0, down, backward}) * g0_k(external_A, {1.0, up, a})
       * g0_k({1.0, up, a}, external_B);
  }

  computed_val = integrand(std::vector<double>{1.0});

  EXPECT_COMPLEX_NEAR(expected_val, computed_val, 1e-17);
}

TEST(integrand_direct, Order_2) { // NOLINT
  model_param_t params;
  g0_model g0{params};
  g0_keldysh_contour_t g0_k{g0};
  auto external_A = gf_index_t{5.0, up, forward};
  auto external_B = gf_index_t{5.0, up, forward};
  integrand_g_direct integrand{g0_k, external_A, external_B};

  dcomplex computed_val;
  dcomplex expected_val;

  //// Order 2 with full determinants (unprecise)
  //gf_index_t A, B;
  //dcomplex det_up, det_down;
  //expected_val = 0;

  //for (keldysh_idx_t a : {forward, backward}) {
  //  for (keldysh_idx_t b : {forward, backward}) {
  //    int sign = (((a + b) % 2) == 0) ? 1 : -1;
  //    A = {1.0, up, a};
  //    B = {2.0, up, b};
  //    det_up =
  //      - (g0_k({1.0, up, forward}, {1.0, up, backward}) * g0_k(external_A, B) - g0_k(external_A, A) * g0_k(A, B)) * g0_k(B, external_B)
  //      + (g0_k(B, A) * g0_k(external_A, B) - g0_k(external_A, A) * g0_k({2.0, up, forward}, {2.0, up, backward})) * g0_k(A, external_B);
  //    A = {1.0, down, a};
  //    B = {2.0, down, b};
  //    det_down = g0_k({1.0, down, forward}, {1.0, up, backward}) * g0_k({2.0, down, forward}, {2.0, up, backward})
  //      - g0_k(B, A) * g0_k(A, B);
  //    expected_val += sign * det_up * det_down;
  //  }
  //}

  //computed_val = integrand(std::vector<double>{1.0, 2.0});

  //EXPECT_COMPLEX_NEAR(expected_val, computed_val, 1e-16); // imag part very close to 0

  // We ignore disconnected diagrams
  expected_val = 0;

  for (keldysh_idx_t a : {forward, backward}) {
    for (keldysh_idx_t b : {forward, backward}) {
      int sign = (((a + b) % 2) == 0) ? 1 : -1;
      expected_val += sign
         * (g0_k({1.0, up, forward}, {1.0, up, backward}) * g0_k(external_A, {2.0, up, b})
               * g0_k({2.0, up, b}, external_B)
            + g0_k(external_A, {1.0, up, a}) * g0_k({2.0, up, forward}, {2.0, up, backward})
               * g0_k({1.0, up, a}, external_B))
         * g0_k({2.0, down, b}, {1.0, down, a}) * g0_k({1.0, down, a}, {2.0, down, b});
      expected_val += sign
         * (g0_k(external_A, {1.0, up, a}) * g0_k({1.0, up, a}, {2.0, up, b}) * g0_k({2.0, up, b}, external_B)
            + g0_k({2.0, up, b}, {1.0, up, a}) * g0_k(external_A, {2.0, up, b}) * g0_k({1.0, up, a}, external_B))
         * (g0_k({1.0, down, forward}, {1.0, down, backward}) * g0_k({2.0, down, forward}, {2.0, down, backward})
            - g0_k({2.0, down, b}, {1.0, down, a}) * g0_k({1.0, down, a}, {2.0, down, b}));
    }
  }

  computed_val = integrand(std::vector<double>{1.0, 2.0});

  EXPECT_COMPLEX_NEAR(expected_val, computed_val, 1e-16); // real part very close to 0
}

MAKE_MAIN; // NOLINT
