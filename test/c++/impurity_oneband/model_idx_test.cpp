#include "keldy/common.hpp"
#include "keldy/impurity_oneband/model.hpp"
#include <gtest/gtest.h>
#include <itertools/itertools.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(contour_pt_test, compare_equality) { // NOLINT
  contour_pt_t a{1.0, forward, 0};
  contour_pt_t b{1.0, forward, 0};
  EXPECT_TRUE(a == b);

  contour_pt_t c{2.0, forward, 0};
  EXPECT_FALSE(a == c);

  contour_pt_t d{1.0, backward, 0};
  EXPECT_FALSE(a == d);

  contour_pt_t e{1.0, forward, 1};
  EXPECT_FALSE(a == e);
}

TEST(contour_pt_test, compare_3way_equality) { // NOLINT
  std::vector<double> times{-1.0, 0.0, 1.0};
  std::vector<int> timeorder{-1, 0, 1};

  for (auto t : times) {
    for (auto to : timeorder) {
      contour_pt_t a{t, forward, to};
      contour_pt_t b{t, forward, to};
      EXPECT_TRUE(compare_3way(a, b) == 0);
      EXPECT_TRUE(compare_3way(b, a) == 0);

      contour_pt_t c{t, backward, to};
      contour_pt_t d{t, backward, to};
      EXPECT_TRUE(compare_3way(c, d) == 0);
      EXPECT_TRUE(compare_3way(d, c) == 0);
    }
  }
}

TEST(contour_pt_test, compare_3way_kidx) { // NOLINT
  std::vector<double> times{-2.0, -1.0, 0.0, 1.0, 2.0};
  std::vector<int> timeorder{-2, -1, 0, 1, 2};

  for (auto [t0, t1] : itertools::product(times, times)) {
    for (auto [to0, to1] : itertools::product(timeorder, timeorder)) {
      contour_pt_t a{t0, forward, to0};
      contour_pt_t b{t1, backward, to1};
      EXPECT_TRUE(compare_3way(a, b) < 0);
      EXPECT_TRUE(compare_3way(b, a) > 0);
    }
  }
}

TEST(contour_pt_test, compare_3way_times) { // NOLINT
  std::vector<int> timeorder{-2, -1, 0, 1, 2};
  for (auto [to0, to1] : itertools::product(timeorder, timeorder)) {
    {
      contour_pt_t a{1.0, forward, to0};
      contour_pt_t b{2.0, forward, to1};
      EXPECT_TRUE(compare_3way(a, b) < 0);
      EXPECT_TRUE(compare_3way(b, a) > 0);
    }
    {
      contour_pt_t a{-1.0, forward, to0};
      contour_pt_t b{-2.0, forward, to1};
      EXPECT_TRUE(compare_3way(a, b) > 0);
      EXPECT_TRUE(compare_3way(b, a) < 0);
    }
    {
      contour_pt_t c{2.0, backward, to0};
      contour_pt_t d{1.0, backward, to1};
      EXPECT_TRUE(compare_3way(c, d) < 0);
      EXPECT_TRUE(compare_3way(d, c) > 0);
    }
    {
      contour_pt_t c{-2.0, backward, to0};
      contour_pt_t d{-1.0, backward, to1};
      EXPECT_TRUE(compare_3way(c, d) > 0);
      EXPECT_TRUE(compare_3way(d, c) < 0);
    }
  }
}

TEST(contour_pt_test, compare_3way_timeorder) { // NOLINT
  std::vector<double> times{-2.0, -1.0, 0.0, 1.0, 2.0};
  for (auto t0 : times) {
    {
      contour_pt_t a{t0, forward, 0};
      contour_pt_t b{t0, forward, 1};
      EXPECT_TRUE(compare_3way(a, b) < 0);
      EXPECT_TRUE(compare_3way(b, a) > 0);
    }
    {
      contour_pt_t a{t0, forward, 0};
      contour_pt_t b{t0, forward, -1};
      EXPECT_TRUE(compare_3way(a, b) > 0);
      EXPECT_TRUE(compare_3way(b, a) < 0);
    }
    {
      contour_pt_t c{t0, backward, 1};
      contour_pt_t d{t0, backward, 0};
      EXPECT_TRUE(compare_3way(c, d) < 0);
      EXPECT_TRUE(compare_3way(d, c) > 0);
    }
    {
      contour_pt_t c{t0, backward, -1};
      contour_pt_t d{t0, backward, 0};
      EXPECT_TRUE(compare_3way(c, d) > 0);
      EXPECT_TRUE(compare_3way(d, c) < 0);
    }  }
}

MAKE_MAIN; // NOLINT
