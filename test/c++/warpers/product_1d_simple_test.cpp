#include "keldy/warpers/product_1d_simple.hpp"
#include "./warper_tests_common.hpp"

#include <triqs/test_tools/arrays.hpp>

using namespace keldy;
using namespace keldy::warpers;
using namespace triqs::arrays;

TEST(SimpleProduct1DWarper, NointerpIdentityConstructor2) { // NOLINT
  double const t_max = 1.5;
  auto const warper = warper_product_1d_simple_t([t_max](double /**/) -> double { return 1. / t_max; },
                                                 [t_max](double u) -> double { return u / t_max; },
                                                 [t_max](double l) -> double { return t_max * l; }, t_max);

  basic_test_warper_at_order_1(warper, t_max, 1e-14);
  basic_test_warper_multidim(warper, t_max, 1e-14);

  std::function<double(double)> l_from_u = [t_max](double ui) { return ui / t_max; };
  function_test_warper(
     warper, t_max, [](std::vector<double> const & /**/) -> double { return 1.; }, vectorize(l_from_u),
     [t_max](std::vector<double> const &li) -> double { return std::pow(t_max, li.size()); });
}

TEST(SimpleProduct1DWarper, NointerpExponentialConstructor) { // NOLINT
  double const t_max = 5.;
  auto const warper = make_product_1d_simple_exponential(t_max, 2.0);

  basic_test_warper_at_order_1(warper, t_max, 1e-9);
  basic_test_warper_multidim(warper, t_max, 1e-9);

  function_test_warper(
     warper, t_max,
     [](std::vector<double> const &ui) -> double {
       return std::exp(-std::accumulate(ui.cbegin(), ui.cend(), 0.) / 2.);
     },
     [t_max](std::vector<double> const &ui) -> std::vector<double> {
       std::vector<double> li = {};
       li.reserve(ui.size());
       for (double u : ui) {
         li.push_back((std::exp(-u / 2.) - 1.) / (std::exp(-t_max / 2.) - 1.));
       }
       return li;
     },
     [t_max](std::vector<double> const &li) -> double {
       double product = 1.;
       for (double l : li) {
         product *= (std::exp(-t_max / 2.) - 1.) * l + 1.;
       }
       return std::pow(2. * (1. - std::exp(-t_max / 2.)), li.size()) / product;
     },
     1e-10, 1e-10, 1e-6);
}

TEST(SimpleProduct1DWarper, NointerpFlatTail) { // NOLINT
  double const t_max = 2.0;
  auto f1 = [](double u) -> double {
    if (u <= 1) {
      return 3. / 2. - u;
    }
    return 0;
  };
  auto f1_integrated = [](double u) -> double {
    if (u < 1) {
      return u * (3 - u) / 2.;
    }
    return 1.;
  };
  auto f1_integrated_inv = [](double l) -> double { return 3. / 2. - std::sqrt(9. / 4. - 2 * l); };

  auto const warper = warper_product_1d_simple_t(f1, f1_integrated, f1_integrated_inv, t_max, false);

  basic_test_warper_at_order_1(warper, t_max, 1e-14);
  basic_test_warper_multidim(warper, t_max, 1e-14);
}

MAKE_MAIN // NOLINT
