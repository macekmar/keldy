#include "keldy/warpers/product_1d_simple.hpp"
#include "./warper_tests_common.hpp"

#include <triqs/test_tools/arrays.hpp>

template <typename T, typename S>
std::vector<T> operator/(std::vector<T> const vector, S const scalar) {
  std::vector<T> output = {};
  for (auto i = vector.begin(); i != vector.end(); ++i) {
    output.push_back(*i / scalar);
  }
  return output;
};

using namespace keldy;
using namespace triqs::arrays;

/// --------------------------------------------------------

TEST(SimpleProduct1DWarper, Default) { // NOLINT
  double const t_max = 1.5;
  auto warper = warper_product_1d_simple_t(t_max);

  basic_test_warper_at_order_1(warper, t_max);
  basic_test_warper_multidim(warper, t_max);

  function_test_warper(
     warper, t_max, [](std::vector<double> const ui) -> double { return 1.; },
     [t_max](std::vector<double> const ui) -> std::vector<double> { return ui / t_max; },
     [t_max](std::vector<double> const li) -> double { return std::pow(t_max, li.size()); });
}

TEST(SimpleProduct1DWarper, IdentityConstructor1) { // NOLINT
  double const t_max = 1.5;
  auto const warper = warper_product_1d_simple_t([](double u) -> double { return 1.; }, t_max, 1e2);

  basic_test_warper_at_order_1(warper, t_max);
  basic_test_warper_multidim(warper, t_max);

  function_test_warper(
     warper, t_max, [](std::vector<double> const ui) -> double { return 1.; },
     [t_max](std::vector<double> const ui) -> std::vector<double> { return ui / t_max; },
     [t_max](std::vector<double> const li) -> double { return std::pow(t_max, li.size()); });
}

TEST(SimpleProduct1DWarper, IdentityConstructor2) { // NOLINT
  double const t_max = 1.5;
  auto const warper =
     warper_product_1d_simple_t([](double u) -> double { return 1.; }, [](double u) -> double { return u; },
                                [](double l) -> double { return l; }, t_max, 1e2);

  basic_test_warper_at_order_1(warper, t_max);
  basic_test_warper_multidim(warper, t_max);

  function_test_warper(
     warper, t_max, [](std::vector<double> const ui) -> double { return 1.; },
     [t_max](std::vector<double> const ui) -> std::vector<double> { return ui / t_max; },
     [t_max](std::vector<double> const li) -> double { return std::pow(t_max, li.size()); });
}

TEST(SimpleProduct1DWarper, ExponentialConstructor1) { // NOLINT
  double const t_max = 5.;
  auto const warper = warper_product_1d_simple_t([](double u) -> double { return std::exp(-u / 2.); }, t_max, 1e5);

  basic_test_warper_at_order_1(warper, t_max, 1e-9);
  basic_test_warper_multidim(warper, t_max, 1e-9);

  function_test_warper(
     warper, t_max,
     [](std::vector<double> const ui) -> double { return std::exp(-std::accumulate(ui.cbegin(), ui.cend(), 0.) / 2.); },
     [t_max](std::vector<double> const ui) -> std::vector<double> {
       std::vector<double> li = {};
       for (double u : ui) {
         li.push_back((std::exp(-u / 2.) - 1.) / (std::exp(-t_max / 2.) - 1.));
       }
       return li;
     },
     [t_max](std::vector<double> const li) -> double {
       double product = 1.;
       for (double l : li) {
         product *= (std::exp(-t_max / 2.) - 1.) * l + 1.;
       }
       return std::pow(2. * (1. - std::exp(-t_max / 2.)), li.size()) / product;
     },
     1e-10, 1e-10, 1e-6);
}

TEST(SimpleProduct1DWarper, ExponentialConstructor2) { // NOLINT
  double const t_max = 5.;
  auto const warper =
     warper_product_1d_simple_t([](double u) -> double { return std::exp(-u / 2.); },
                                [](double u) -> double { return 2. * (1 - std::exp(-u / 2.)); },
                                [](double l) -> double { return -2. * std::log(1 - l / 2.); }, t_max, 1e5);

  basic_test_warper_at_order_1(warper, t_max, 1e-9);
  basic_test_warper_multidim(warper, t_max, 1e-9);

  function_test_warper(
     warper, t_max,
     [](std::vector<double> const ui) -> double { return std::exp(-std::accumulate(ui.cbegin(), ui.cend(), 0.) / 2.); },
     [t_max](std::vector<double> const ui) -> std::vector<double> {
       std::vector<double> li = {};
       for (double u : ui) {
         li.push_back((std::exp(-u / 2.) - 1.) / (std::exp(-t_max / 2.) - 1.));
       }
       return li;
     },
     [t_max](std::vector<double> const li) -> double {
       double product = 1.;
       for (double l : li) {
         product *= (std::exp(-t_max / 2.) - 1.) * l + 1.;
       }
       return std::pow(2. * (1. - std::exp(-t_max / 2.)), li.size()) / product;
     },
     1e-10, 1e-10, 1e-6);
}

MAKE_MAIN; // NOLINT
