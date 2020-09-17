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

/// Check if any element of `container` equals `value`
template <typename T, typename V>
void expect_contain(T const &container, V value) {
  bool found = false;

  for (auto const &v : container) {
    if (v == value) {
      found = true;
      break;
    }
  }

  EXPECT_TRUE(found);
}

TEST(g0_model, MiddleTimePoint1) { // NOLINT
  model_param_t params;
  params.nr_time_points_gf = 1000;
  g0_model g0{g0_model_omega{params}, true};

  expect_contain(g0.g0_lesser[up].mesh(), 0.);
  expect_contain(g0.g0_greater[up].mesh(), 0.);
}

TEST(g0_model, MiddleTimePoint2) { // NOLINT
  model_param_t params;
  params.nr_time_points_gf = 1001;
  g0_model g0{g0_model_omega{params}, true};

  expect_contain(g0.g0_lesser[up].mesh(), 0.);
  expect_contain(g0.g0_greater[up].mesh(), 0.);
}

TEST(g0_keldysh_adaptor, Initialize) { // NOLINT
  model_param_t params;
  g0_model g0{g0_model_omega{params}, true};
  g0_keldysh_contour_t g0_k{g0};
}

TEST(g0_keldysh_adaptor, LesserGreater) { // NOLINT
  model_param_t params;
  params.eps_d = 0.0;
  params.alpha = 0.0;
  g0_model g0{g0_model_omega{params}, false};
  g0_keldysh_contour_t g0_k{g0};

  double const tol = 1e-12;

  // time, spin, k_idx, timesplit_n, orbital
  gf_index_t const A = {2.0, up, forward, 0, 0};
  gf_index_t const B = {2.0, up, forward, 1, 0};
  gf_index_t const C = {3.0, up, forward, 0, 0};
  gf_index_t const D = {3.0, up, backward, 0, 0};
  gf_index_t const E = {2.0, up, backward, 1, 0};
  gf_index_t const F = {2.0, up, backward, 0, 0};

  EXPECT_COMPLEX_NEAR(g0_k(A, B), g0.g0_lesser[up](A.contour.time - B.contour.time)(0, 0), tol);
  EXPECT_COMPLEX_NEAR(g0_k(B, C), g0.g0_lesser[up](B.contour.time - C.contour.time)(0, 0), tol);
  EXPECT_COMPLEX_NEAR(g0_k(C, D), g0.g0_lesser[up](C.contour.time - D.contour.time)(0, 0), tol);
  EXPECT_COMPLEX_NEAR(g0_k(D, E), g0.g0_lesser[up](D.contour.time - E.contour.time)(0, 0), tol);
  EXPECT_COMPLEX_NEAR(g0_k(E, F), g0.g0_lesser[up](E.contour.time - F.contour.time)(0, 0), tol);

  EXPECT_COMPLEX_NEAR(g0_k(F, E), g0.g0_greater[up](F.contour.time - E.contour.time)(0, 0), tol);
  EXPECT_COMPLEX_NEAR(g0_k(E, D), g0.g0_greater[up](E.contour.time - D.contour.time)(0, 0), tol);
  EXPECT_COMPLEX_NEAR(g0_k(D, C), g0.g0_greater[up](D.contour.time - C.contour.time)(0, 0), tol);
  EXPECT_COMPLEX_NEAR(g0_k(C, B), g0.g0_greater[up](C.contour.time - B.contour.time)(0, 0), tol);
  EXPECT_COMPLEX_NEAR(g0_k(B, A), g0.g0_greater[up](B.contour.time - A.contour.time)(0, 0), tol);
}

TEST(g0_keldysh_adaptor, Orbital) { // NOLINT
  model_param_t params;
  params.eps_d = 0.0;
  params.alpha = 0.0;
  g0_model g0{g0_model_omega{params}, true};
  g0_keldysh_contour_t g0_k{g0};

  double const tol = 1e-12;

  // time, spin, k_idx, timesplit_n, orbital
  // dot:
  gf_index_t const A = {-10.0, up, forward, 0, 0};
  gf_index_t const B = {-1.0, up, backward, 0, 0};
  // lead:
  gf_index_t const C = {1.0, up, forward, 0, 1};
  gf_index_t const D = {10.0, up, backward, 0, 1};

  EXPECT_COMPLEX_NEAR(g0_k(A, C), g0.g0_lesser[up](A.contour.time - C.contour.time)(0, 1), tol);
  EXPECT_COMPLEX_NEAR(g0_k(A, D), g0.g0_lesser[up](A.contour.time - D.contour.time)(0, 1), tol);
  EXPECT_COMPLEX_NEAR(g0_k(B, C), g0.g0_greater[up](B.contour.time - C.contour.time)(0, 1), tol);
  EXPECT_COMPLEX_NEAR(g0_k(B, D), g0.g0_greater[up](B.contour.time - D.contour.time)(0, 1), tol);

  EXPECT_COMPLEX_NEAR(g0_k(C, A), g0.g0_greater[up](C.contour.time - A.contour.time)(1, 0), tol);
  EXPECT_COMPLEX_NEAR(g0_k(D, A), g0.g0_greater[up](D.contour.time - A.contour.time)(1, 0), tol);
  EXPECT_COMPLEX_NEAR(g0_k(C, B), g0.g0_lesser[up](C.contour.time - B.contour.time)(1, 0), tol);
  EXPECT_COMPLEX_NEAR(g0_k(D, B), g0.g0_lesser[up](D.contour.time - B.contour.time)(1, 0), tol);

  // In anticipation of larger impurities
  EXPECT_COMPLEX_NEAR(g0_k(C, D), g0.g0_lesser[up](C.contour.time - D.contour.time)(1, 1), tol);
  EXPECT_COMPLEX_NEAR(g0_k(D, C), g0.g0_greater[up](D.contour.time - C.contour.time)(1, 1), tol);
}

TEST(g0_keldysh_adaptor, WithAlpha) { // NOLINT
  using namespace std::complex_literals;

  model_param_t params;
  params.eps_d = 0.0;
  params.alpha = 0.8;
  g0_model g0{g0_model_omega{params}, true};
  g0_keldysh_contour_t g0_k{g0};

  double const tol = 1e-12;

  // time, spin, k_idx, timesplit_n, orbital
  gf_index_t const A = {2.0, up, forward, 0, 0};
  gf_index_t const B = {2.0, up, forward, 1, 0};
  gf_index_t const C = {2.0, up, forward, 0, 1};

  EXPECT_COMPLEX_NEAR(g0_k(A, A), g0.g0_lesser[up](0.0)(0, 0) - 1.0i * 0.8, tol);
  EXPECT_COMPLEX_NEAR(g0_k(A, A, true), g0.g0_lesser[up](0.0)(0, 0) - 1.0i * 0.8, tol);
  EXPECT_COMPLEX_NEAR(g0_k(B, B), g0.g0_lesser[up](0.0)(0, 0) - 1.0i * 0.8, tol);
  EXPECT_COMPLEX_NEAR(g0_k(C, C), g0.g0_lesser[up](0.0)(1, 1) - 1.0i * 0.8, tol);
}

TEST(g0_keldysh_adaptor, WithoutAlpha) { // NOLINT
  model_param_t params;
  params.eps_d = 0.0;
  params.alpha = 0.8;
  g0_model g0{g0_model_omega{params}, true};
  g0_keldysh_contour_t g0_k{g0};

  double const tol = 1e-12;

  // time, spin, k_idx, timesplit_n, orbital
  gf_index_t A = {2.0, up, forward, 0, 0};
  gf_index_t B = {3.0, up, forward, 0, 0};
  EXPECT_COMPLEX_NEAR(g0_k(A, B), g0.g0_lesser[up](-1.0)(0, 0), tol);

  A = {2.0, up, forward, 0, 0};
  B = {2.0, down, forward, 0, 0};
  EXPECT_COMPLEX_NEAR(g0_k(A, B), 0., tol);

  A = {2.0, up, forward, 0, 0};
  B = {2.0, up, backward, 0, 0};
  EXPECT_COMPLEX_NEAR(g0_k(A, B), g0.g0_lesser[up](0.0)(0, 0), tol);

  A = {2.0, up, forward, 0, 0};
  B = {2.0, up, forward, 1, 0};
  EXPECT_COMPLEX_NEAR(g0_k(A, B), g0.g0_lesser[up](0.0)(0, 0), tol);

  A = {2.0, up, forward, 0, 0};
  B = {2.0, up, forward, 0, 1};
  EXPECT_COMPLEX_NEAR(g0_k(A, B), g0.g0_lesser[up](0.0)(0, 1), tol);
}

MAKE_MAIN // NOLINT
