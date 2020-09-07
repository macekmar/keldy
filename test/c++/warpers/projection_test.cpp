#include "keldy/warpers/plasma_uv.hpp"
#include "keldy/warpers/projection.hpp"
#include "./warper_tests_common.hpp"

#include <triqs/test_tools/arrays.hpp>

using namespace keldy;
using namespace keldy::warpers;

using namespace triqs::arrays;

/// --------------------------------------------------------

TEST(ProjectionWarper, OneToOne) { // NOLINT
  auto integrand = [](std::vector<double> const xi) -> std::vector<dcomplex> {
    double output = 1.;
    for (double x : xi) {
      output *= std::sin(3 * x);
    }
    return {output, 1};
  };

  auto proj_warper = warper_projection_t(integrand, 4, int(1e3), int(1e3), 0.06, false);

  // TODO: that is not very good precision
  basic_test_warper_at_order_1(proj_warper, 1.0, 1e-6);
  basic_test_warper_multidim(proj_warper, 1.0, 1e-6);
}

TEST(ProjectionWarper, OptimizeSigma) { // NOLINT
  double const tmax = 1.;
  double const tol = 1e-6;
  auto integrand = [](std::vector<double> const xi) -> std::vector<dcomplex> {
    double output = 1.;
    for (double x : xi) {
      output *= std::sin(3 * x);
    }
    return {output, 1};
  };

  // TODO: get rid of this warper in tests?
  auto cst = []([[maybe_unused]] double x) -> double { return 1.; };
  warper_product_1d_simple_interp_nearest_t warper = {cst, tmax, 10};

  auto warped_integrand = [&warper, &integrand, tmax](std::vector<double> const li) -> std::vector<dcomplex> {
    auto [vi, jac] = warper.map_forward(li);
    auto ui = ui_from_vi(tmax, vi);
    return {integrand(ui)[0] * jac, 1};
  };

  auto proj_warper = warper_projection_t(warped_integrand, 1, int(1e2), int(1e3), 0.1, true);

  /// TODO: fix minimizer and adjust reference value:
  //EXPECT_NEAR(proj_warper.get_sigmas()[0], 0.06754688329, tol);

  auto proj_warper_higher_order = warper_projection_t(warped_integrand, 3, int(1e2), int(1e3), 0.1, true);
}

TEST(ProjectionWarper, Values) { // NOLINT
  double const tmax = 1.;
  double const tol = 5e-6;
  auto integrand = [](std::vector<double> const xi) -> std::vector<dcomplex> {
    double output = 1.;
    for (double x : xi) {
      output *= std::sin(3 * x);
    }
    return {output, 1};
  };

  // TODO: get rid of this warper in tests?
  auto cst = []([[maybe_unused]] double x) -> double { return 1.; };
  warper_product_1d_simple_interp_nearest_t warper = {cst, tmax, 10};

  auto warped_integrand = [&warper, &integrand, tmax](std::vector<double> const li) -> std::vector<dcomplex> {
    auto [vi, jac] = warper.map_forward(li);
    auto ui = ui_from_vi(tmax, vi);
    return {integrand(ui)[0] * jac, 1};
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
/// TODO: add more substantial tests

MAKE_MAIN // NOLINT
