#pragma once

#include <triqs/test_tools/arrays.hpp>
#include <vector>

/*
 * Some functions to help test iterables, like std::vector or std::array
 */

template <typename T>
size_t get_iterable_size(T const &a) {
  size_t size = 0;
  for (auto i = a.cbegin(); i < a.cend(); ++i) {
    size++;
  }
  return size;
}

/* Check element-wise strict equality
 * Use like:
 * EXPECT_TRUE(are_iterable_eq(a, b));
 */
template <typename T>
bool are_iterable_eq(T const &a, T const &b) {
  auto size_a = get_iterable_size(a);
  auto size_b = get_iterable_size(b);
  if (size_a == size_b) {
    auto i = a.cbegin();
    auto j = b.cbegin();
    for (; i < a.cend(); ++i, ++j) {
      if (*i != *j) {
        std::cout << "All elements are not equal, e.g.: " << *i << " != " << *j << " (diff = " << *i - *j << ")"
                  << std::endl;
        return false;
      }
    }
  } else {
    std::cout << "Iterables have different sizes: " << size_a << " != " << size_b << std::endl;
    return false;
  }
  return true;
}

/* Check element-wise closeness
 * Use like:
 * EXPECT_TRUE(are_iterable_near(a, b));
 */
template <typename T>
bool are_iterable_near(T const &a, T const &b, double abs_tol) {
  auto size_a = get_iterable_size(a);
  auto size_b = get_iterable_size(b);
  if (size_a == size_b) {
    auto i = a.cbegin();
    auto j = b.cbegin();
    for (; i < a.cend(); ++i, ++j) {
      if (std::abs(*i - *j) >= abs_tol) {
        std::cout << "All elements are not close, e.g.: " << *i << " - " << *j << " = " << *i - *j << " >= " << abs_tol
                  << std::endl;
        return false;
      }
    }
  } else {
    std::cout << "Iterables have different sizes: " << size_a << " != " << size_b << std::endl;
    return false;
  }
  return true;
}

template <typename T>
bool are_iterable_near(T const &a, T const &b, std::vector<double> const &abs_tol) {
  auto size_a = get_iterable_size(a);
  auto size_b = get_iterable_size(b);
  auto size_c = get_iterable_size(abs_tol);
  if (size_a == size_b) {
    if (size_c != size_a) {
      std::cout << "Vector abs_tol error is not size compatible with a / b" << std::endl;
      return false;
    }
    auto i = a.cbegin();
    auto j = b.cbegin();
    auto k = abs_tol.cbegin();
    for (; i < a.cend(); ++i, ++j, ++k) {
      if (std::abs(*i - *j) >= *k) {
        std::cout << "All elements are not close, e.g.: " << *i << " - " << *j << " = " << *i - *j << " >= " << *k
                  << std::endl;
        return false;
      }
    }
  } else {
    std::cout << "Iterables have different sizes: " << size_a << " != " << size_b << std::endl;
    return false;
  }
  return true;
}

/* Check element-wize almost equality
 * TODO: in case of failure, line number refers to this file, which is pointless
 */
template <typename T>
void expect_iterable_double_eq(T const &a, T const &b) {
  auto size_a = get_iterable_size(a);
  auto size_b = get_iterable_size(b);
  EXPECT_EQ(size_a, size_b);
  if (size_a == size_b) {
    auto i = a.cbegin();
    auto j = b.cbegin();
    for (; i < a.cend(); ++i, ++j) {
      EXPECT_DOUBLE_EQ(*i, *j);
    }
  }
}
