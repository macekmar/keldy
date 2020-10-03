#include "keldy/warpers/product_1d_simple.hpp"
#include "./warper_tests_common.hpp"

#include <triqs/test_tools/arrays.hpp>

using namespace keldy;
using namespace keldy::warpers;

using namespace triqs::arrays;

/// --------------------------------------------------------

TEST(SimpleProduct1DWarper, InterpNearest_Default) { // NOLINT
  double const t_max = 1.5;
  auto warper = warper_product_1d_simple_interp_nearest_t([](double /**/) { return 1.0; }, t_max, 8);

  basic_test_warper_at_order_1(warper, t_max, 1e-10, true);
  basic_test_warper_multidim(warper, t_max, 1e-10, true);

  std::function<double(double)> l_from_u = [t_max](double const &ui) { return (t_max - ui) / t_max; };
  function_test_warper(
     warper, t_max, [](std::vector<double> const /**/) -> double { return 1.; }, vectorize(l_from_u),
     [t_max](std::vector<double> const &li) -> double { return std::pow(t_max, li.size()); });
}

TEST(SimpleProduct1DWarper, InterpNearest_IdentityConstructor) { // NOLINT
  double const t_max = 1.5;
  auto const warper = warper_product_1d_simple_interp_nearest_t([](double /**/) -> double { return 1.; }, t_max, 1.0e2);

  basic_test_warper_at_order_1(warper, t_max, 1e-10, true);
  basic_test_warper_multidim(warper, t_max, 1e-10, true);

  std::function<double(double)> l_from_u = [t_max](double const &ui) { return (t_max - ui) / t_max; };
  function_test_warper(
     warper, t_max, [](std::vector<double> const /**/) -> double { return 1.; }, vectorize(l_from_u),
     [t_max](std::vector<double> const &li) -> double { return std::pow(t_max, li.size()); });
}

TEST(SimpleProduct1DWarper, InterpNearest_ExponentialConstructor) { // NOLINT
  double const t_max = 5.;
  auto const warper =
     warper_product_1d_simple_interp_nearest_t([](double u) -> double { return std::exp(-u / 2.); }, t_max, int(1e3));

  basic_test_warper_at_order_1(warper, t_max, 1e-9, true);
  basic_test_warper_multidim(warper, t_max, 1e-9, true);

  std::function<double(double)> l_from_u = [t_max](double const &ui) {
    return (std::exp(-t_max / 2.) - std::exp(-ui / 2.)) / std::expm1(-t_max / 2.);
  };

  function_test_warper(
     warper, t_max,
     [](std::vector<double> const &ui) -> double {
       return std::exp(-std::accumulate(ui.cbegin(), ui.cend(), 0.) / 2.);
     },
     vectorize(l_from_u),
     [t_max](std::vector<double> const &li) -> double {
       double product = 1.;
       for (double l : li) {
         product *= -std::expm1(-t_max / 2.) * l + std::exp(-t_max / 2.);
       }
       return std::pow(-2. * std::expm1(-t_max / 2.), li.size()) / product;
     },
     1e-3, 3e-6, 1.1e-3, true);
  // errors are ~ 1/delta, 1/delta^2, 1/delta
}
TEST(SimpleProduct1DWarper, InterpNearest_ExponentialConstructor_LongTail) { // NOLINT
  double const t_max = 500.;
  auto const warper =
     warper_product_1d_simple_interp_nearest_t([](double u) -> double { return std::exp(-u / 2.); }, t_max, int(1e3));

  basic_test_warper_at_order_1(warper, t_max, 1e-9, true);
  basic_test_warper_multidim(warper, t_max, 1e-9, true);
}

/// --------------------------------------------------------

TEST(SimpleProduct1DWarper, InterpHybrid_Default) { // NOLINT
  double const t_max = 1.5;
  auto warper = warper_product_1d_simple_interp_hybrid_t([](double /**/) { return 1.0; }, t_max, 8);

  basic_test_warper_at_order_1(warper, t_max, 1e-10, true);
  basic_test_warper_multidim(warper, t_max, 1e-10, true);

  std::function<double(double)> l_from_u = [t_max](double const &ui) { return (t_max - ui) / t_max; };
  function_test_warper(
     warper, t_max, [](std::vector<double> const /**/) -> double { return 1.; }, vectorize(l_from_u),
     [t_max](std::vector<double> const &li) -> double { return std::pow(t_max, li.size()); });
}

TEST(SimpleProduct1DWarper, InterpHybrid_IdentityConstructor) { // NOLINT
  double const t_max = 1.5;
  auto const warper = warper_product_1d_simple_interp_hybrid_t([](double /**/) -> double { return 1.; }, t_max, 1.0e2);

  basic_test_warper_at_order_1(warper, t_max, 1e-10, true);
  basic_test_warper_multidim(warper, t_max, 1e-10, true);

  std::function<double(double)> l_from_u = [t_max](double const &ui) { return (t_max - ui) / t_max; };
  function_test_warper(
     warper, t_max, [](std::vector<double> const /**/) -> double { return 1.; }, vectorize(l_from_u),
     [t_max](std::vector<double> const &li) -> double { return std::pow(t_max, li.size()); });
}

TEST(SimpleProduct1DWarper, InterpHybrid_ExponentialConstructor) { // NOLINT
  double const t_max = 5.;
  auto const warper = warper_product_1d_simple_interp_hybrid_t([](double u) -> double { return std::exp(-u / 2.); },
                                                               t_max, 5 * int(1e3));

  basic_test_warper_at_order_1(warper, t_max, 1e-9, true);
  basic_test_warper_multidim(warper, t_max, 1e-9, true);

  std::function<double(double)> l_from_u = [t_max](double const &ui) {
    return (std::exp(-t_max / 2.) - std::exp(-ui / 2.)) / std::expm1(-t_max / 2.);
  };

  function_test_warper(
     warper, t_max,
     [](std::vector<double> const &ui) -> double {
       return std::exp(-std::accumulate(ui.cbegin(), ui.cend(), 0.) / 2.);
     },
     vectorize(l_from_u),
     [t_max](std::vector<double> const &li) -> double {
       double product = 1.;
       for (double l : li) {
         product *= -std::expm1(-t_max / 2.) * l + std::exp(-t_max / 2.);
       }
       return std::pow(-2. * std::expm1(-t_max / 2.), li.size()) / product;
     },
     1e-10, 1e-10, 1e-6);
}

TEST(SimpleProduct1DWarper, InterpHybrid_ExponentialConstructor_LongTail) { // NOLINT
  double const t_max = 400.;
  auto const warper =
     warper_product_1d_simple_interp_hybrid_t([](double u) -> double { return std::exp(-u / 2.); }, t_max, int(1e4));

  basic_test_warper_at_order_1(warper, t_max, 1e-6, true);
  basic_test_warper_multidim(warper, t_max, 1e-6, true);
}

MAKE_MAIN // NOLINT
