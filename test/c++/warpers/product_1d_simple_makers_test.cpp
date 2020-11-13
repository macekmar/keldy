#include "keldy/warpers/product_1d_simple.hpp"
#include "./warper_tests_common.hpp"

#include <triqs/test_tools/arrays.hpp>

using namespace keldy::warpers;

TEST(SimpleProduct1DWarperMaker, Exponential) { // NOLINT
  auto warper_1 = make_product_1d_simple_exponential(10.0, 1.5);
  auto warper_2 = make_product_1d_simple_exponential(1.0, 1.8);
  auto warper_3 = make_product_1d_simple_exponential(20.0, 1.1);

  double const tol = 1e-9;
  basic_test_warper_at_order_1(warper_1, 10.0, tol);
  basic_test_warper_at_order_1(warper_2, 1.0, tol);
  basic_test_warper_at_order_1(warper_3, 20.0, tol);

  basic_test_warper_multidim(warper_1, 10.0, tol);
  basic_test_warper_multidim(warper_2, 1.0, tol);
  basic_test_warper_multidim(warper_3, 20.0, tol);

  EXPECT_DOUBLE_EQ(warper_1.jacobian_forward({2., 4.}) / warper_1.jacobian_forward({3., 4.}),
                   std::exp((3. - 2.) / 1.5));
}

TEST(SimpleProduct1DWarperMaker, Inverse) { // NOLINT
  auto warper_1 = make_product_1d_simple_inverse(10.0, 1.5);
  auto warper_2 = make_product_1d_simple_inverse(1.0, 1.8);
  auto warper_3 = make_product_1d_simple_inverse(20.0, 1.1);

  double const tol = 1e-13;
  basic_test_warper_at_order_1(warper_1, 10.0, tol);
  basic_test_warper_at_order_1(warper_2, 1.0, tol);
  basic_test_warper_at_order_1(warper_3, 20.0, tol);

  basic_test_warper_multidim(warper_1, 10.0, tol);
  basic_test_warper_multidim(warper_2, 1.0, tol);
  basic_test_warper_multidim(warper_3, 20.0, tol);

  EXPECT_DOUBLE_EQ(warper_1.jacobian_forward({2., 4.}) / warper_1.jacobian_forward({3., 4.}), (3. + 1.5) / (2. + 1.5));
}

TEST(SimpleProduct1DWarperMaker, InverseSquare) { // NOLINT
  auto warper_1 = make_product_1d_simple_inverse_square(10.0, 1.5);
  auto warper_2 = make_product_1d_simple_inverse_square(1.0, 1.8);
  auto warper_3 = make_product_1d_simple_inverse_square(20.0, 1.1);

  double const tol = 1e-13;
  basic_test_warper_at_order_1(warper_1, 10.0, tol);
  basic_test_warper_at_order_1(warper_2, 1.0, tol);
  basic_test_warper_at_order_1(warper_3, 20.0, tol);

  basic_test_warper_multidim(warper_1, 10.0, tol);
  basic_test_warper_multidim(warper_2, 1.0, tol);
  basic_test_warper_multidim(warper_3, 20.0, tol);

  EXPECT_DOUBLE_EQ(warper_1.jacobian_forward({2., 4.}) / warper_1.jacobian_forward({3., 4.}),
                   std::pow((3. + 1.5) / (2. + 1.5), 2));
}

TEST(SimpleProduct1DWarperMaker, InversePower) { // NOLINT
  auto warper_1 = make_product_1d_simple_inverse_power(2, 10.0, 1.5);
  auto warper_2 = make_product_1d_simple_inverse_power(3, 1.0, 1.8);
  auto warper_3 = make_product_1d_simple_inverse_power(4, 20.0, 1.1);

  double const tol = 1e-10;
  basic_test_warper_at_order_1(warper_1, 10.0, tol);
  basic_test_warper_at_order_1(warper_2, 1.0, tol);
  basic_test_warper_at_order_1(warper_3, 20.0, tol);

  basic_test_warper_multidim(warper_1, 10.0, tol);
  basic_test_warper_multidim(warper_2, 1.0, tol);
  basic_test_warper_multidim(warper_3, 20.0, tol);

  EXPECT_DOUBLE_EQ(warper_1.jacobian_forward({2., 4.}) / warper_1.jacobian_forward({3., 4.}),
                   std::pow((3. + 1.5) / (2. + 1.5), 2));
  EXPECT_DOUBLE_EQ(warper_2.jacobian_forward({0.2, 0.4}) / warper_2.jacobian_forward({0.3, 0.4}),
                   std::pow((0.3 + 1.8) / (0.2 + 1.8), 3));
  EXPECT_DOUBLE_EQ(warper_3.jacobian_forward({2., 4.}) / warper_3.jacobian_forward({3., 4.}),
                   std::pow((3. + 1.1) / (2. + 1.1), 4));
}

MAKE_MAIN // NOLINT
