#include "keldy/warpers/product_1d_simple.hpp"
#include "./warper_tests_common.hpp"

#include <triqs/test_tools/arrays.hpp>

using namespace keldy::warpers;

TEST(SimpleProduct1DWarperMaker, Exponential) { // NOLINT
  auto warper_1 = make_product_1d_simple_exponential_nointerp(10.0, 1.5);
  auto warper_2 = make_product_1d_simple_exponential_nointerp(1.0, 1.8);
  auto warper_3 = make_product_1d_simple_exponential_nointerp(20.0, 1.1);
}


TEST(SimpleProduct1DWarperMaker, Inverse) { // NOLINT
  auto warper_1 = make_product_1d_simple_inverse_nointerp(10.0, 1.5);
  auto warper_2 = make_product_1d_simple_inverse_nointerp(1.0, 1.8);
  auto warper_3 = make_product_1d_simple_inverse_nointerp(20.0, 1.1);
}


TEST(SimpleProduct1DWarperMaker, InverseSquare) { // NOLINT
  auto warper_1 = make_product_1d_simple_inverse_square_nointerp(10.0, 1.5);
  auto warper_2 = make_product_1d_simple_inverse_square_nointerp(1.0, 1.8);
  auto warper_3 = make_product_1d_simple_inverse_square_nointerp(20.0, 1.1);
}



MAKE_MAIN; // NOLINT
