#include "keldy/warpers/plasma_uv.hpp"
#include "keldy/warpers/projection.hpp"
#include "keldy/warpers/make_warper_from_proj.hpp"
#include "./warper_tests_common.hpp"

#include <triqs/test_tools/arrays.hpp>

using namespace keldy;
using namespace keldy::warpers;

using namespace triqs::arrays;

/// --------------------------------------------------------

TEST(ProjectionWarper, OneToOne) { // NOLINT
  auto integrand = [](std::vector<double> const xi) -> double {
    double output = 1.;
    for (double x : xi) {
      output *= std::sin(3 * x);
    }
    return output;
  };

  auto proj_warper = warper_projection_t(integrand, 4, int(1e3), int(1e3), 0.06, false);

  // TODO: that is not very good precision
  basic_test_warper_at_order_1(proj_warper, 1.0, 1e-6);
  basic_test_warper_multidim(proj_warper, 1.0, 1e-6);
}

TEST(ProjectionWarper, OneToOne_Maker) { // NOLINT
  auto integrand = [](std::vector<double> const xi) -> double {
    double output = 1.;
    for (double x : xi) {
      output *= std::sin(3 * x);
    }
    return output;
  };

  auto proj_warper = make_warper_from_proj(integrand, 4, int(1e3), int(1e3), 0.06);

  // TODO: that is not very good precision
  basic_test_warper_at_order_1(proj_warper, 1.0, 1e-6, true);
  basic_test_warper_multidim(proj_warper, 1.0, 1e-6, true);
}

TEST(ProjectionWarper, OptimizeSigma) { // NOLINT
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
  auto cst = []([[maybe_unused]] double x) -> double { return 1.; };
  warper_product_1d_simple_interp_nearest_t warper = {cst, tmax, 10};

  auto warped_integrand = [&warper, &integrand, tmax](std::vector<double> const li) -> double {
    auto [vi, jac] = warper.map_forward(li);
    auto ui = ui_from_vi(tmax, vi);
    return integrand(ui) * jac;
  };

  auto proj_warper = warper_projection_t(warped_integrand, 1, int(1e2), int(1e3), 0.1, true);

  /// TODO: fix minimizer and adjust reference value:
  //EXPECT_NEAR(proj_warper.get_sigmas()[0], 0.06754688329, tol);

  auto proj_warper_higher_order = warper_projection_t(warped_integrand, 3, int(1e2), int(1e3), 0.1, true);
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
  auto cst = []([[maybe_unused]] double x) -> double { return 1.; };
  warper_product_1d_simple_interp_nearest_t warper = {cst, tmax, 10};

  auto warped_integrand = [&warper, &integrand, tmax](std::vector<double> const li) -> double {
    auto [vi, jac] = warper.map_forward(li);
    auto ui = ui_from_vi(tmax, vi);
    return integrand(ui) * jac;
  };

  auto proj_warper = warper_projection_t(warped_integrand, 1, int(1e2), int(1e2), 0.06754688329, false);

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

  EXPECT_NEAR(proj_warper.jacobian_forward({0.3}), 0.841834, tol);
  EXPECT_NEAR(proj_warper.jacobian_forward({0.5}), 0.977436, tol);
  EXPECT_NEAR(proj_warper.jacobian_forward({0.9}), 0.237891, tol);
  EXPECT_NEAR(proj_warper.ui_from_li({0.3})[0], 0.341224, tol);
  EXPECT_NEAR(proj_warper.ui_from_li({0.5})[0], 0.47522, tol);
  EXPECT_NEAR(proj_warper.ui_from_li({0.9})[0], 0.77635, tol);
  EXPECT_NEAR(proj_warper.li_from_ui({0.3})[0], 0.243646, tol);
  EXPECT_NEAR(proj_warper.li_from_ui({0.5})[0], 0.538089, tol);
  EXPECT_NEAR(proj_warper.li_from_ui({0.9})[0], 0.98085, tol);
  EXPECT_NEAR(proj_warper.jacobian_reverse({0.3}), 0.709755, tol);
  EXPECT_NEAR(proj_warper.jacobian_reverse({0.5}), 0.650334, tol);
  EXPECT_NEAR(proj_warper.jacobian_reverse({0.9}), 1.08937, tol);
}

TEST(ProjectionWarper, Values_Maker) { // NOLINT
  double const tmax = 1.;
  double const tol = 1e-3;
  auto integrand = [](std::vector<double> const xi) -> double {
    double output = 1.;
    for (double x : xi) {
      output *= std::sin(3 * x);
    }
    return output;
  };

  // TODO: get rid of this warper in tests?
  auto cst = []([[maybe_unused]] double x) -> double { return 1.; };
  warper_product_1d_simple_interp_nearest_t warper = {cst, tmax, 10};

  auto warped_integrand = [&warper, &integrand, tmax](std::vector<double> const li) -> double {
    auto [vi, jac] = warper.map_forward(li);
    auto ui = ui_from_vi(tmax, vi);
    return integrand(ui) * jac;
  };

  auto proj_warper1 = warper_projection_t(warped_integrand, 1, int(5e3), int(1e5), 0.06754688329, false);
  auto proj_warper2 = make_warper_from_proj(warped_integrand, 1, int(5e3), int(1e5), 0.06754688329);

  EXPECT_NEAR(proj_warper1.ui_from_li({0.3})[0], proj_warper2.ui_from_li({1. - 0.3})[0], tol);
  EXPECT_NEAR(proj_warper1.ui_from_li({0.5})[0], proj_warper2.ui_from_li({1. - 0.5})[0], tol);
  EXPECT_NEAR(proj_warper1.ui_from_li({0.9})[0], proj_warper2.ui_from_li({1. - 0.9})[0], tol);
  EXPECT_NEAR(proj_warper1.li_from_ui({0.3})[0], 1. - proj_warper2.li_from_ui({0.3})[0], tol);
  EXPECT_NEAR(proj_warper1.li_from_ui({0.5})[0], 1. - proj_warper2.li_from_ui({0.5})[0], tol);
  EXPECT_NEAR(proj_warper1.li_from_ui({0.9})[0], 1. - proj_warper2.li_from_ui({0.9})[0], tol);
  EXPECT_NEAR(proj_warper1.jacobian_reverse({0.3}), proj_warper2.jacobian_reverse({1. - 0.3}), tol);
  EXPECT_NEAR(proj_warper1.jacobian_reverse({0.5}), proj_warper2.jacobian_reverse({1. - 0.5}), tol);
  EXPECT_NEAR(proj_warper1.jacobian_reverse({0.9}), proj_warper2.jacobian_reverse({1. - 0.9}), tol);
}
/// TODO: add more substantial tests

MAKE_MAIN // NOLINT
