#pragma once

#include "../tests_std_containers.hpp"

#include <gtest/gtest.h>
#include <functional>

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
// with 0 mapped to 0 and t_max mapped to 1.
// The uv warper will fail at this test.

using warper_func_t = std::function<double(std::vector<double>)>;
using warper_map_t = std::function<std::vector<double>(std::vector<double>)>;

template <typename W>
inline void basic_test_warper_at_order_1(W const &warper, double const t_max, double const accuracy) {
  EXPECT_GT(t_max, 0.);

  /// boundaries are almost equal
  EXPECT_NEAR(warper.li_from_ui({0.})[0], 0., accuracy);
  EXPECT_NEAR(warper.li_from_ui({t_max})[0], 1., accuracy);
  EXPECT_NEAR(warper.ui_from_li({0.})[0], 0., accuracy);
  EXPECT_NEAR(warper.ui_from_li({1.})[0], t_max, accuracy);

  /// boundaries are respected
  for (double u : {0., 0.01 * t_max, 0.5 * t_max, 0.99 * t_max, t_max}) {
    EXPECT_GE(warper.li_from_ui({u})[0], -accuracy);
    EXPECT_LE(warper.li_from_ui({u})[0] - 1., accuracy);
  }
  for (double l : {0., 0.01, 0.5, 0.99, 1.}) {
    EXPECT_GE(warper.ui_from_li({l})[0], -accuracy);
    EXPECT_LE(warper.ui_from_li({l})[0] - t_max, accuracy);
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
inline void basic_test_warper_multidim(W const &warper, double const t_max, double const accuracy) {
  EXPECT_GT(t_max, 0.);
  std::vector<double> ui = {};
  std::vector<double> li = {};
  std::vector<double> li_plus = {};

  /// boundaries are almost equal
  EXPECT_TRUE(are_iterable_near(warper.li_from_ui({0., 0., 0., 0.}), {0., 0., 0., 0.}, accuracy));
  EXPECT_TRUE(are_iterable_near(warper.ui_from_li({0., 0., 0., 0.}), {0., 0., 0., 0.}, accuracy));
  EXPECT_TRUE(are_iterable_near(warper.li_from_ui({t_max, t_max, t_max, t_max}), {1., 1., 1., 1.}, accuracy));
  EXPECT_TRUE(are_iterable_near(warper.ui_from_li({1., 1., 1., 1.}), {t_max, t_max, t_max, t_max}, accuracy));

  /// boundaries are respected
  ui = {0.7 * t_max, t_max, 0.3 * t_max, 0.};
  EXPECT_TRUE(vector_is_bounded(warper.li_from_ui(ui), -accuracy, 1. + accuracy));

  /// difference in l is resolved in u (i.e. non-flatness of the l(u) curve)
  for (double l : {0.0001, 0.01, 0.5, 0.99, 0.9999}) {
    li = {0.7, l, 0.3, 0.5};
    li_plus = {0.7, l + 1e-9, 0.3, 0.5};
    EXPECT_NE(warper.ui_from_li(li), warper.ui_from_li(li_plus));
    EXPECT_NE(warper.map_reverse(li).first, warper.map_reverse(li_plus).first);
  }

  /// reverse mapping
  for (double u = 0.1 * t_max; u < t_max; u += 0.1 * t_max) {
    ui = {0.7 * t_max, t_max, u, 0.};
    EXPECT_TRUE(are_iterable_near(warper.ui_from_li(warper.li_from_ui(ui)), ui, accuracy));
    EXPECT_TRUE(are_iterable_near(warper.map_reverse(warper.map_forward(ui).first).first, ui, accuracy));

  }
  for (double l = 0.1; l < 1.; l += 0.1) {
    li = {0.7, 1., l, 0.};
    EXPECT_TRUE(are_iterable_near(warper.li_from_ui(warper.ui_from_li(li)), li, accuracy));
    EXPECT_TRUE(are_iterable_near(warper.map_forward(warper.map_reverse(li).first).first, li, accuracy));
  }
};

/// check that the different methods of the warper correspond to the functions provided
template <typename W>
inline void function_test_warper(W const &warper, double const t_max, warper_func_t const &f,
                                 warper_map_t const &li_from_ui_ref, warper_func_t const &jacobian_ref,
                                 double accuracy_f = 1e-10, double accuracy_map = 1e-10, double accuracy_jac = 1e-10,
                                 bool use_accuracy_relative = false) {
  std::vector<double> li = {};

  auto do_test = [&](std::vector<double> const ui) -> void {
    std::vector<double> t_max_vec(ui.size(), t_max);

    li = li_from_ui_ref(ui);

    std::vector<double> accuracy_map_vec(li.size(), accuracy_map);

    if (use_accuracy_relative) {
      accuracy_f *= warper.jacobian_forward(ui) / warper.jacobian_forward(t_max_vec);
      for (int i = 0; i < li.size(); i++) {
        if (li.at(i) != 0.0) {
          accuracy_map_vec.at(i) *= li.at(i);
        }
      }
      accuracy_jac *= warper.jacobian_reverse(li);
    }

    EXPECT_NEAR(warper.jacobian_forward(ui) / warper.jacobian_forward(t_max_vec), f(ui) / f(t_max_vec), accuracy_f);
    EXPECT_TRUE(are_iterable_near(warper.li_from_ui(ui), li_from_ui_ref(ui), accuracy_map_vec));

    EXPECT_NEAR(warper.jacobian_reverse(li), jacobian_ref(li), accuracy_jac);
    // EXPECT_NEAR(warper.jacobian_forward(li), 1.0 / jacobian_ref(li), accuracy_jac);
    EXPECT_GT(warper.jacobian_reverse(li), 0.); // has to be != 0.

    // Repeat test with map forward / backwards
    auto [li_from_ui_vec, jac_for] = warper.map_forward(ui);
    EXPECT_NEAR(jac_for / warper.jacobian_forward(t_max_vec), f(ui) / f(t_max_vec), accuracy_f);
    EXPECT_TRUE(are_iterable_near(li_from_ui_vec, li_from_ui_ref(ui), accuracy_map_vec));

    auto [ui_from_li_vec, jac_rev] = warper.map_reverse(li);
    EXPECT_NEAR(jac_rev, jacobian_ref(li), accuracy_jac);
    EXPECT_GT(jac_rev, 0.); // has to be != 0.
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
