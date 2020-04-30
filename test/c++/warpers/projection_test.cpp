#include "keldy/warpers/plasma_uv.hpp"
#include "keldy/warpers/plasma_projection.hpp"
#include "./warper_tests_common.hpp"

#include <triqs/test_tools/arrays.hpp>

using namespace keldy;
using namespace keldy::warpers;

using namespace triqs::arrays;

/// --------------------------------------------------------

// FIX THIS TEST

TEST(ProjectionWarper, OneToOne) { // NOLINT
  auto integrand = [](std::vector<double> const xi) -> double {
    double output = 1.;
    for (double x : xi) {
      output *= std::sin(3 * x);
    }
    return output;
  };

  auto cst = [](double x) -> double { return 1.; };
  warper_product_1d_t warper = {{cst, cst, cst, cst}, 1., 10};

  EXPECT_EQ(0,1); // fix rest of test

  // auto proj_warper = warper_plasma_projection_t(integrand, warper, 4, 1e3);

  // basic_test_warper_at_order_1(proj_warper, 1.);
  // basic_test_warper_multidim(proj_warper, 1.);
}

/// TODO: add more substantial tests

MAKE_MAIN; // NOLINT
