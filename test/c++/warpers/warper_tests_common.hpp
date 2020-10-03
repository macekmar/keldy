#pragma once

#include "../tests_std_containers.hpp"

#include <gtest/gtest.h>
#include <functional>

template <typename T>
std::function<std::vector<T>(std::vector<T>)> vectorize(std::function<T(T)> const &f) {
  auto new_func = [&f](std::vector<T> x) -> std::vector<T> {
    std::vector<T> y = {};
    y.reserve(x.size());
    for (auto i = x.cbegin(); i != x.cend(); ++i) {
      y.push_back(f(*i));
    }
    return y;
  };
  return new_func;
}

template <typename T>
bool vector_is_bounded(std::vector<T> const &x, T const x_min, T const x_max) {
  bool output = true;
  for (auto const &elt : x) {
    output = output && (x_min <= elt) && (elt <= x_max);
  }
  return output;
}

///--------------------------- Warper tests -----------------------------
// 'basic' tests are for warpers of type u -> l,
// They check that ... TODO
// The uv warper will fail at this test.

using warper_func_t = std::function<double(std::vector<double>)>;
using warper_map_t = std::function<std::vector<double>(std::vector<double>)>;

template <typename W>
inline void basic_test_warper_at_order_1(W const &warper, double const u_max, double const accuracy,
                                         bool const negative = false, double const l_min = 0, double const l_max = 1) {
  double const u_min = 0;
  ASSERT_LT(u_min, u_max);
  ASSERT_LT(l_min, l_max);

  std::vector<double> list{1e-9, 0.0001, 0.001, 0.1, 0.2, 0.5, 0.7, 0.999, 0.9999, 1. - 1e-9};
  auto u_list = list; // copy
  auto l_list = list; // copy
  auto stretch_u = [u_min, u_max](double x) -> double { return u_min + x * (u_max - u_min); };
  auto stretch_l = [l_min, l_max](double x) -> double { return l_min + x * (l_max - l_min); };
  std::transform(u_list.begin(), u_list.end(), u_list.begin(), stretch_u);
  std::transform(l_list.begin(), l_list.end(), l_list.begin(), stretch_l);

  // if warper is decreasing, monotony has to be checked in reverse
  if (negative) {
    std::reverse(u_list.begin(), u_list.end());
    std::reverse(l_list.begin(), l_list.end());
  }

  /// u -> l boundaries
  if (not negative) {
    EXPECT_NEAR(warper.li_from_ui({u_min})[0], l_min, accuracy);
    EXPECT_NEAR(warper.li_from_ui({u_max})[0], l_max, accuracy);
  } else {
    EXPECT_NEAR(warper.li_from_ui({u_min})[0], l_max, accuracy);
    EXPECT_NEAR(warper.li_from_ui({u_max})[0], l_min, accuracy);
  }
  for (auto const u : u_list) {
    double const l = warper.li_from_ui({u})[0];
    EXPECT_GE(l, l_min - accuracy);
    EXPECT_LE(l, l_max + accuracy);
  }

  /// u -> l monotony
  // cannot check slope, because continuity is only piecewise (rolle thm does not apply)
  double l_tmp_max = -std::numeric_limits<double>::max();
  for (auto const u : u_list) {
    double const l = warper.li_from_ui({u})[0];
    EXPECT_GE(l, l_tmp_max);
    l_tmp_max = l;
    /// very close u values
    if (u + 1e-9 <= u_max) {
      if (not negative) {
        EXPECT_LE(l, warper.li_from_ui({u + 1e-9})[0]);
      } else {
        EXPECT_GE(l, warper.li_from_ui({u + 1e-9})[0]);
      }
    }
  }

  /// l -> u boundaries
  if (not negative) {
    EXPECT_NEAR(warper.ui_from_li({l_min})[0], u_min, accuracy); // only one side
  } else {
    EXPECT_NEAR(warper.ui_from_li({l_max})[0], u_min, accuracy); // only one side
  }
  for (auto const l : l_list) {
    double const u = warper.ui_from_li({l})[0];
    EXPECT_GE(u, u_min - accuracy);
    EXPECT_LE(u, u_max + accuracy);
  }

  /// l -> u strict monotony
  double u_tmp_max = -std::numeric_limits<double>::max();
  for (auto const l : l_list) {
    double const u = warper.ui_from_li({l})[0];
    EXPECT_GT(u, u_tmp_max); // strict
    u_tmp_max = u;
    /// very close l values
    if (l + 1e-9 <= l_max) {
      if (not negative) {
        EXPECT_LT(u, warper.ui_from_li({l + 1e-9})[0]); // strict
      } else {
        EXPECT_GT(u, warper.ui_from_li({l + 1e-9})[0]); // strict
      }
    }
  }

  /// inversion u -> l -> u
  bool tested_once = false;
  for (auto const u : u_list) {
    if (warper.jacobian_forward({u}) > 0.) { // if model function is not zero
      tested_once = true;
      EXPECT_NEAR(u, warper.ui_from_li(warper.li_from_ui({u}))[0], accuracy);
    }
  }
  EXPECT_TRUE(tested_once);

  /// inversion l -> u -> l
  for (auto const l : l_list) {
    EXPECT_NEAR(l, warper.li_from_ui(warper.ui_from_li({l}))[0], accuracy);
  }

  /// map reverse
  for (auto const l : l_list) {
    auto [ui, jac] = warper.map_reverse({l});
    EXPECT_EQ(ui[0], warper.ui_from_li({l})[0]);
    EXPECT_EQ(jac, warper.jacobian_reverse({l}));
  }

  /// map forward
  for (auto const u : u_list) {
    auto [li, jac] = warper.map_forward({u});
    EXPECT_EQ(li[0], warper.li_from_ui({u})[0]);
    EXPECT_EQ(jac, warper.jacobian_forward({u}));
  }
}

template <typename W>
inline void basic_test_warper_multidim(W const &warper, double const u_max, double const accuracy,
                                       bool const negative = false, double const l_min = 0, double const l_max = 1) {
  double const u_min = 0;
  ASSERT_LT(u_min, u_max);
  ASSERT_LT(l_min, l_max);

  std::vector<double> ui = {};
  std::vector<double> li = {};
  std::vector<double> li_plus = {};

  std::vector<double> list{1e-9, 0.0001, 0.001, 0.1, 0.2, 0.5, 0.7, 0.999, 0.9999, 1. - 1e-9};
  auto u_list = list; // copy
  auto l_list = list; // copy
  auto stretch_u = [u_min, u_max](double x) -> double { return u_min + x * (u_max - u_min); };
  auto stretch_l = [l_min, l_max](double x) -> double { return l_min + x * (l_max - l_min); };
  std::transform(u_list.begin(), u_list.end(), u_list.begin(), stretch_u);
  std::transform(l_list.begin(), l_list.end(), l_list.begin(), stretch_l);

  /// boundaries are almost equal (except l = 1 ---> u <= t_max)
  std::vector<double> vec_u_min = {u_min, u_min, u_min, u_min};
  std::vector<double> vec_u_max = {u_max, u_max, u_max, u_max};
  std::vector<double> vec_l_min = {l_min, l_min, l_min, l_min};
  std::vector<double> vec_l_max = {l_max, l_max, l_max, l_max};
  if (not negative) {
    EXPECT_TRUE(are_iterable_near(warper.li_from_ui(vec_u_min), vec_l_min, accuracy));
    EXPECT_TRUE(are_iterable_near(warper.li_from_ui(vec_u_max), vec_l_max, accuracy));
    EXPECT_TRUE(are_iterable_near(warper.ui_from_li(vec_l_min), vec_u_min, accuracy));
  } else {
    EXPECT_TRUE(are_iterable_near(warper.li_from_ui(vec_u_min), vec_l_max, accuracy));
    EXPECT_TRUE(are_iterable_near(warper.li_from_ui(vec_u_max), vec_l_min, accuracy));
    EXPECT_TRUE(are_iterable_near(warper.ui_from_li(vec_l_max), vec_u_min, accuracy));
  }

  /// u -> l boundaries
  for (auto const u : u_list) {
    ui = {stretch_u(0.7), stretch_u(1.), u, stretch_u(0.)};
    EXPECT_PRED3(vector_is_bounded<double>, warper.li_from_ui(ui), l_min - accuracy, l_max + accuracy);
  }

  /// l -> u boundaries
  for (auto const l : l_list) {
    li = {stretch_l(0.3), l, stretch_l(0.9), stretch_l(0.1)};
    EXPECT_PRED3(vector_is_bounded<double>, warper.ui_from_li(li), u_min - accuracy, u_max + accuracy);
  }

  /// different l gives different u (opposite not always true)
  for (auto const l : l_list) {
    if (l + 1e-9 <= l_max) {
      li = {stretch_l(0.7), l, stretch_l(0.3), stretch_l(0.5)};
      li_plus = li;
      li_plus[1] += 1e-9;
      EXPECT_NE(warper.ui_from_li(li), warper.ui_from_li(li_plus));
    }
  }
  // TODO: what of the monotony in multidim warper?

  /// l -> u -> l is identity
  for (auto const l : l_list) {
    li = {stretch_l(0.7), stretch_l(1.), l, stretch_l(0.)};
    EXPECT_TRUE(are_iterable_near(warper.li_from_ui(warper.ui_from_li(li)), li, accuracy));
  }

  /// map reverse
  for (auto const l : l_list) {
    li = {stretch_l(0.7), stretch_l(1.), l, stretch_l(0.)};
    auto [ui, jac] = warper.map_reverse(li);
    EXPECT_PRED2(are_iterable_eq<std::vector<double>>, ui, warper.ui_from_li(li));
    EXPECT_EQ(jac, warper.jacobian_reverse(li));
  }

  /// map forward
  for (auto const u : u_list) {
    ui = {stretch_u(0.7), stretch_u(1.), u, stretch_u(0.)};
    auto [li, jac] = warper.map_forward(ui);
    EXPECT_PRED2(are_iterable_eq<std::vector<double>>, li, warper.li_from_ui(ui));
    EXPECT_EQ(jac, warper.jacobian_forward(ui));
  }
}

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
      accuracy_f *= std::abs(warper.jacobian_forward(ui) / warper.jacobian_forward(t_max_vec));
      for (int i = 0; i < li.size(); i++) {
        if (li.at(i) != 0.0) {
          accuracy_map_vec.at(i) *= std::abs(li.at(i));
        }
      }
      accuracy_jac *= std::abs(warper.jacobian_reverse(li));
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
}
