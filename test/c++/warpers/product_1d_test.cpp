#include "keldy/warpers/product_1d.hpp"
#include "./warper_tests_common.hpp"

#include <triqs/test_tools/arrays.hpp>

using namespace keldy;
using namespace keldy::warpers;

using namespace triqs::arrays;

/// --------------------------------------------------------

TEST(Product1DWarper, Default) { // NOLINT
  double const t_max = 1.5;

  std::function<double(double)> const_function = []([[maybe_unused]] double t) { return 1.0; };
  auto warper = warper_product_1d_interp_nearest_t{};
  warper.emplace_back({const_function, t_max, 8});

  /// only order 1 exists
  basic_test_warper_at_order_1(warper, t_max, 1e-14, true);
}

TEST(Product1DWarper, Identity) { // NOLINT
  double const t_max = 1.5;
  auto cst = []([[maybe_unused]] double u) -> double { return 1.; };

  auto warper = warper_product_1d_interp_nearest_t{};
  warper.emplace_back({cst, t_max, int(1e2)});
  warper.emplace_back({cst, t_max, int(1e2)});
  warper.emplace_back({cst, t_max, int(1e2)});
  warper.emplace_back({cst, t_max, int(1e2)});

  basic_test_warper_at_order_1(warper, t_max, 1e-14, true);
  basic_test_warper_multidim(warper, t_max, 1e-14, true);

  std::function<double(double)> l_from_u = [t_max](double const &ui) { return (t_max - ui) / t_max; };
  function_test_warper(
     warper, t_max, []([[maybe_unused]] std::vector<double> const ui) -> double { return 1.; }, vectorize(l_from_u),
     [t_max](std::vector<double> const li) -> double { return std::pow(t_max, li.size()); });
}

TEST(Product1DWarper, AlternateInverseAndCube) { // NOLINT
  double const t_max = 5.;
  auto f1 = [](double u) -> double { return 1. / (2. + u); };
  auto f2 = [](double u) -> double { return 1. / ((3. + u) * (3. + u) * (3. + u)); };

  auto warper = warper_product_1d_interp_nearest_t{};
  warper.emplace_back({f1, t_max, int(1e5)});
  warper.emplace_back({f2, t_max, int(1e5)});
  warper.emplace_back({f1, t_max, int(1e5)});
  warper.emplace_back({f2, t_max, int(1e5)});

  basic_test_warper_at_order_1(warper, t_max, 1e-9, true);
  basic_test_warper_multidim(warper, t_max, 1e-9, true);

  auto f = [&f1, &f2](std::vector<double> const ui) -> double {
    double output = 1.;
    for (size_t i = 0; i < ui.size(); ++i) {
      output *= (i % 2 == 0) ? f1(ui[i]) : f2(ui[i]);
    }
    return output;
  };

  double const norm1 = std::log(1 + t_max / 2.);
  double const norm2 = 1. / ((3. + t_max) * (3. + t_max)) - 1. / 9.;

  auto li_from_ui = [norm1, norm2](std::vector<double> const ui) {
    std::vector<double> output = {};
    for (size_t i = 0; i < ui.size(); ++i) {
      if (i % 2 == 0) {
        output.push_back((norm1 - std::log(1 + ui[i] / 2.)) / norm1);
      } else {
        output.push_back((norm2 - 1. / ((3. + ui[i]) * (3. + ui[i])) + 1. / 9.) / norm2);
      }
    }
    return output;
  };

  auto jac = [norm1, norm2](std::vector<double> li) -> double {
    double output = 1.;
    for (size_t i = 0; i < li.size(); ++i) {
      if (i % 2 == 0) {
        output *= 2. * norm1 * std::exp((1 - li[i]) * norm1);
      } else {
        output *= -norm2 / std::pow((1 - li[i]) * norm2 + 1. / 9., 1.5) / 2.;
      }
    }
    return output;
  };

  function_test_warper(warper, t_max, f, li_from_ui, jac, 1e-5, 1e-10, 1e-5, true);
  // errors are ~ 1/delta, 1/delta^2, 1/delta
}

TEST(Product1DWarper, AlternateInverseAndCubeHybrid) { // NOLINT
  double const t_max = 5.;
  auto f1 = [](double u) -> double { return 1. / (2. + u); };
  auto f2 = [](double u) -> double { return 1. / ((3. + u) * (3. + u) * (3. + u)); };

  auto warper = warper_product_1d_interp_hybrid_t{};
  warper.emplace_back({f1, t_max, int(1e5)});
  warper.emplace_back({f2, t_max, int(1e5)});
  warper.emplace_back({f1, t_max, int(1e5)});
  warper.emplace_back({f2, t_max, int(1e5)});

  basic_test_warper_at_order_1(warper, t_max, 1e-9, true);
  basic_test_warper_multidim(warper, t_max, 1e-9, true);

  auto f = [&f1, &f2](std::vector<double> const ui) -> double {
    double output = 1.;
    for (size_t i = 0; i < ui.size(); ++i) {
      output *= (i % 2 == 0) ? f1(ui[i]) : f2(ui[i]);
    }
    return output;
  };

  double const norm1 = std::log(1 + t_max / 2.);
  double const norm2 = 1. / ((3. + t_max) * (3. + t_max)) - 1. / 9.;

  auto li_from_ui = [norm1, norm2](std::vector<double> const ui) {
    std::vector<double> output = {};
    for (size_t i = 0; i < ui.size(); ++i) {
      if (i % 2 == 0) {
        output.push_back((norm1 - std::log(1 + ui[i] / 2.)) / norm1);
      } else {
        output.push_back((norm2 - 1. / ((3. + ui[i]) * (3. + ui[i])) + 1. / 9.) / norm2);
      }
    }
    return output;
  };

  auto jac = [norm1, norm2](std::vector<double> li) -> double {
    double output = 1.;
    for (size_t i = 0; i < li.size(); ++i) {
      if (i % 2 == 0) {
        output *= 2. * norm1 * std::exp((1 - li[i]) * norm1);
      } else {
        output *= -norm2 / std::pow((1 - li[i]) * norm2 + 1. / 9., 1.5) / 2.;
      }
    }
    return output;
  };

  function_test_warper(warper, t_max, f, li_from_ui, jac, 1e-14, 1e-11, 1e-11, true);
  // errors are ~ 1/delta, 1/delta^2, 1/delta
}

MAKE_MAIN // NOLINT
