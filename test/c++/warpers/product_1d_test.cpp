#include "keldy/warpers/product_1d.hpp"
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

TEST(Product1DWarper, Default) { // NOLINT
  double const t_max = 1.5;
  auto warper = warper_product_1d_t(t_max);

  /// only order 1 exists
  basic_test_warper_at_order_1(warper, t_max);
}

TEST(Product1DWarper, Identity) { // NOLINT
  double const t_max = 1.5;
  auto cst = [](double u) -> double { return 1.; };
  auto const warper = warper_product_1d_t({cst, cst, cst, cst}, t_max, 1e2);

  basic_test_warper_at_order_1(warper, t_max);
  basic_test_warper_multidim(warper, t_max);

  function_test_warper(
     warper, t_max, [](std::vector<double> const ui) -> double { return 1.; },
     [t_max](std::vector<double> const ui) -> std::vector<double> { return ui / t_max; },
     [t_max](std::vector<double> const li) -> double { return std::pow(t_max, li.size()); });
}

TEST(Product1DWarper, AlternateInverseAndCube) { // NOLINT
  double const t_max = 5.;
  auto f1 = [](double u) -> double { return 1. / (2. + u); };
  auto f2 = [](double u) -> double { return 1. / ((3. + u) * (3. + u) * (3. + u)); };
  auto const warper = warper_product_1d_t({f1, f2, f1, f2}, t_max, 1e5);

  basic_test_warper_at_order_1(warper, t_max, 1e-9);
  basic_test_warper_multidim(warper, t_max, 1e-9);

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
        output.push_back(std::log(1 + ui[i] / 2.) / norm1);
      } else {
        output.push_back((1. / ((3. + ui[i]) * (3. + ui[i])) - 1. / 9.) / norm2);
      }
    }
    return output;
  };

  auto jac = [norm1, norm2](std::vector<double> li) -> double {
    double output = 1.;
    for (size_t i = 0; i < li.size(); ++i) {
      if (i % 2 == 0) {
        output *= 2. * norm1 * std::exp(li[i] * norm1);
      } else {
        output *= -norm2 / std::pow(li[i] * norm2 + 1. / 9., 1.5) / 2.;
      }
    }
    return output;
  };

  function_test_warper(warper, t_max, f, li_from_ui, jac, 1e-10, 1e-10, 1e-6);
}

MAKE_MAIN; // NOLINT
