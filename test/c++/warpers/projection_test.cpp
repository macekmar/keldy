#include "keldy/warpers/plasma_uv.hpp"
#include "keldy/warpers/projection.hpp"
#include "./warper_tests_common.hpp"

#include <triqs/test_tools/arrays.hpp>

using namespace keldy;
using namespace keldy::warpers;

using namespace triqs::arrays;

/// --------------------------------------------------------

TEST(ProjectionWarper, OneToOne) { // NOLINT
  double const tmax = 1.;
  auto integrand = [](std::vector<double> const xi) -> double {
    double output = 1.;
    for (double x : xi) {
      output *= std::sin(3 * x);
    }
    return output;
  };

  auto proj_warper = warper_projection_t(integrand, tmax, 4, int(1e3), int(1e3), 0.06, false);

  // TODO: that is not very good precision
  basic_test_warper_at_order_1(proj_warper, tmax, 1e-6);
  basic_test_warper_multidim(proj_warper, tmax, 1e-6);
}

TEST(ProjectionWarper, SigmaValue) { // NOLINT
  double const tmax = 1.;
  double const tol = 1e-6;
  auto integrand = [](std::vector<double> const xi) -> double {
    double output = 1.;
    for (double x : xi) {
      output *= std::sin(3 * x);
    }
    return output;
  };

  // TODO: get rid of this warper in tests?
  auto cst = [](double x) -> double { return 1.; };
  warper_product_1d_simple_t warper = {cst, tmax, 10};

  auto warped_integrand = [&warper, &integrand, tmax](std::vector<double> const li) -> double {
    auto [vi, jac] = warper.map_forward(li);
    auto ui = ui_from_vi(tmax, vi);
    return integrand(ui) * jac;
  };

  auto proj_warper = warper_projection_t(warped_integrand, tmax, 1, int(1e2), int(1e2), 0.1, true);

  EXPECT_NEAR(proj_warper.get_sigmas()[0], 0.06754688329, tol);
}

TEST(ProjectionWarper, Values) { // NOLINT
  double const tmax = 1.;
  double const tol = 5e-6;
  auto integrand = [](std::vector<double> const xi) -> double {
    double output = 1.;
    for (double x : xi) {
      output *= std::sin(3 * x);
    }
    return output;
  };

  // TODO: get rid of this warper in tests?
  auto cst = [](double x) -> double { return 1.; };
  warper_product_1d_simple_t warper = {cst, tmax, 10};

  auto warped_integrand = [&warper, &integrand, tmax](std::vector<double> const li) -> double {
    auto [vi, jac] = warper.map_forward(li);
    auto ui = ui_from_vi(tmax, vi);
    return integrand(ui) * jac;
  };

  auto proj_warper = warper_projection_t(warped_integrand, tmax, 1, int(1e2), int(1e2), 0.06754688329, false);

  //std::cout << proj_warper.jacobian_forward({0.3}) << std::endl;
  //std::cout << proj_warper.jacobian_forward({0.5}) << std::endl;
  //std::cout << proj_warper.jacobian_forward({0.9}) << std::endl;
  //std::cout << proj_warper.ui_from_li({0.3})[0] << std::endl;
  //std::cout << proj_warper.ui_from_li({0.5})[0] << std::endl;
  //std::cout << proj_warper.ui_from_li({0.9})[0] << std::endl;
  //std::cout << proj_warper.li_from_ui({0.3})[0] << std::endl;
  //std::cout << proj_warper.li_from_ui({0.5})[0] << std::endl;
  //std::cout << proj_warper.li_from_ui({0.9})[0] << std::endl;
  //std::cout << proj_warper.jacobian_reverse({0.3}) << std::endl;
  //std::cout << proj_warper.jacobian_reverse({0.5}) << std::endl;
  //std::cout << proj_warper.jacobian_reverse({0.9}) << std::endl;

  EXPECT_NEAR(proj_warper.jacobian_forward({0.3}), 0.764648, tol);
  EXPECT_NEAR(proj_warper.jacobian_forward({0.5}), 0.857939, tol);
  EXPECT_NEAR(proj_warper.jacobian_forward({0.9}), 0.273381, tol);
  EXPECT_NEAR(proj_warper.ui_from_li({0.3})[0], 0.341149, tol);
  EXPECT_NEAR(proj_warper.ui_from_li({0.5})[0], 0.481958, tol);
  EXPECT_NEAR(proj_warper.ui_from_li({0.9})[0], 0.790271, tol);
  EXPECT_NEAR(proj_warper.li_from_ui({0.3})[0], 0.242981, tol);
  EXPECT_NEAR(proj_warper.li_from_ui({0.5})[0], 0.527598, tol);
  EXPECT_NEAR(proj_warper.li_from_ui({0.9})[0], 0.974175, tol);
  EXPECT_NEAR(proj_warper.jacobian_reverse({0.3}), 0.720125, tol);
  EXPECT_NEAR(proj_warper.jacobian_reverse({0.5}), 0.661651, tol);
  EXPECT_NEAR(proj_warper.jacobian_reverse({0.9}), 1.21275, tol);
}
/// TODO: add more substantial tests

MAKE_MAIN; // NOLINT
