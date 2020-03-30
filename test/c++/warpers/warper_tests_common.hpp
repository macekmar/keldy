#pragma once

#include <triqs/test_tools/arrays.hpp>

#include <vector>
#include <functional>

///------------ std::vector checks -------------

void expect_vectors_are_equal(std::vector<double> const &x, std::vector<double> const &y) {
  if (x.size() == y.size()) {
    for (auto const &[xi, yi] : itertools::zip(x, y)) {
      EXPECT_DOUBLE_EQ(xi, yi);
    }
  } else {
    std::cout << "Vector sizes do not match: " << x.size() << " != " << y.size() << std::endl;
    ADD_FAILURE();
  }
};

template <typename T>
bool vectors_are_near(std::vector<T> const &x, std::vector<T> const &y, T const abs_err) {
  if (x.size() == y.size()) {
    bool output = true;
    for (auto const &[xi, yi] : itertools::zip(x, y)) {
      output = output && (std::abs(xi - yi) <= abs_err);
    }
    return output;
  } else {
    std::cout << "Vector sizes do not match: " << x.size() << " != " << y.size() << std::endl;
    return false;
  }
};

template <typename T>
bool vector_is_bounded(std::vector<T> const &x, T const x_min, T const x_max) {
  bool output = true;
  for (auto const &elt : x) {
    output = output && (x_min <= elt) && (elt <= x_max);
  }
  return output;
};

///--------------------------- Warper tests -----------------------------
// 'basic' tests are for warpers of type u-> l, they check that it forms a
// one-to-one map between [0, t_max]^n and [0, 1]^n

using warper_func_t = std::function<double(std::vector<double>)>;
using warper_map_t = std::function<std::vector<double>(std::vector<double>)>;

template <typename W>
inline void basic_test_warper_at_order_1(W const &warper, double const t_max, double const accuracy = 1e-10) {
  EXPECT_GT(t_max, 0.);

  /// boundaries are almost equal
  EXPECT_DOUBLE_EQ(warper.li_from_ui({0.})[0], 0.);
  EXPECT_DOUBLE_EQ(warper.li_from_ui({t_max})[0], 1.);
  EXPECT_DOUBLE_EQ(warper.ui_from_li({0.})[0], 0.);
  EXPECT_DOUBLE_EQ(warper.ui_from_li({1.})[0], t_max);

  /// boundaries are respected
  for (double u : {0., 0.01 * t_max, 0.5 * t_max, 0.99 * t_max, t_max}) {
    EXPECT_GE(warper.li_from_ui({u})[0], 0.);
    EXPECT_LE(warper.li_from_ui({u})[0], 1.);
  }
  for (double l : {0., 0.01, 0.5, 0.99, 1.}) {
    EXPECT_GE(warper.ui_from_li({l})[0], 0.);
    EXPECT_LE(warper.ui_from_li({l})[0], t_max);
  }

  /// difference in l is resolved in u (i.e. non-flatness of the l(u) curve)
  for (double l : {0.0001, 0.01, 0.5, 0.99, 0.9999}) {
    // assuming we need 10^9 samples, each having a different coordinate in l:
    EXPECT_NE(warper.ui_from_li({l})[0], warper.ui_from_li({l + 1e-9})[0]);
  }

  /// reverse mapping
  for (double u = 0.1 * t_max; u < t_max; u += 0.1 * t_max) {
    EXPECT_NEAR(warper.ui_from_li(warper.li_from_ui({u}))[0], u, accuracy);
  }
  for (double l = 0.1; l < 1.; l += 0.1) {
    EXPECT_NEAR(warper.li_from_ui(warper.ui_from_li({l}))[0], l, accuracy);
  }
};

template <typename W>
inline void basic_test_warper_multidim(W const &warper, double const t_max, double const accuracy = 1e-10) {
  EXPECT_GT(t_max, 0.);
  std::vector<double> ui = {};
  std::vector<double> li, li_plus = {};

  /// boundaries are almost equal
  expect_vectors_are_equal(warper.li_from_ui({0., 0., 0., 0.}), {0., 0., 0., 0.});
  expect_vectors_are_equal(warper.ui_from_li({0., 0., 0., 0.}), {0., 0., 0., 0.});
  expect_vectors_are_equal(warper.li_from_ui({t_max, t_max, t_max, t_max}), {1., 1., 1., 1.});
  expect_vectors_are_equal(warper.ui_from_li({1., 1., 1., 1.}), {t_max, t_max, t_max, t_max});

  /// boundaries are respected
  ui = {0.7 * t_max, t_max, 0.3 * t_max, 0.};
  EXPECT_TRUE(vector_is_bounded(warper.li_from_ui(ui), 0., 1.));

  /// difference in l is resolved in u (i.e. non-flatness of the l(u) curve)
  for (double l : {0.0001, 0.01, 0.5, 0.99, 0.9999}) {
    li = {0.7, l, 0.3, 0.5};
    li_plus = {0.7, l + 1e-9, 0.3, 0.5};
    EXPECT_NE(warper.ui_from_li(li), warper.ui_from_li(li_plus));
  }

  /// reverse mapping
  for (double u = 0.1 * t_max; u < t_max; u += 0.1 * t_max) {
    ui = {0.7 * t_max, t_max, u, 0.};
    EXPECT_TRUE(vectors_are_near(warper.ui_from_li(warper.li_from_ui(ui)), ui, accuracy));
  }
  for (double l = 0.1; l < 1.; l += 0.1) {
    li = {0.7, 1., l, 0.};
    EXPECT_TRUE(vectors_are_near(warper.li_from_ui(warper.ui_from_li(li)), li, accuracy));
  }
};

/// check that the different methods of the warper correspond to the functions provided
template <typename W>
inline void function_test_warper(W const &warper, double const t_max, warper_func_t const &f,
                                 warper_map_t const &li_from_ui_ref, warper_func_t const &jacobian_ref,
                                 double const accuracy_f = 1e-10, double const accuracy_map = 1e-10,
                                 double const accuracy_jac = 1e-10) {
  std::vector<double> li = {};

  auto do_test = [&](std::vector<double> const ui) -> void {
    EXPECT_NEAR(warper(ui), f(ui), accuracy_f);
    EXPECT_TRUE(vectors_are_near(warper.li_from_ui(ui), li_from_ui_ref(ui), accuracy_map));

    li = li_from_ui_ref(ui);
    EXPECT_NEAR(warper.jacobian(li), jacobian_ref(li), accuracy_jac);
    EXPECT_GT(warper.jacobian(li), 0.); // has to be != 0.
  };

  for (double u = 0.1 * t_max; u < t_max; u += 0.1 * t_max) {
    do_test({u});
    do_test({u, 0.6 * t_max});
    do_test({u, t_max});
    do_test({0.2 * t_max, u, 0.99 * t_max});
    do_test({u, 0., 0.3 * t_max});
    do_test({0.7 * t_max, t_max, u, 0.});
  }
};
