#include "keldy/warpers/plasma_uv.hpp"
#include <keldy/warpers/warpers.hpp>

#include <keldy/common.hpp>

#include <itertools/itertools.hpp>
#include <triqs/gfs.hpp>
#include <triqs/test_tools/arrays.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy;
using namespace keldy::warpers;
using namespace triqs::arrays;

TEST(ViUiMaps, Simple) { // NOLINT
  double t_max = 1.0;
  std::vector<double> ui_times = {0.8, 0.6, 0.2, 0.1};
  std::vector<double> vi_out = vi_from_ui(t_max, ui_times);

  std::vector<double> vi_result = {0.2, 0.2, 0.4, 0.1};

  for (auto const &[v1, v2] : itertools::zip(vi_out, vi_result)) {
    EXPECT_DOUBLE_EQ(v1, v2);
  }
  EXPECT_EQ(ui_from_vi(t_max, vi_out), ui_times);
}

TEST(ViUiMaps, OutOfOrder) { // NOLINT
  double t_max = 10.0;
  std::vector<double> ui_times = {0.0, 2.0, 6.0, 1.0};
  std::vector<double> vi_out = vi_from_ui(t_max, ui_times);

  std::vector<double> vi_result = {4.0, 4.0, 1.0, 1.0};

  for (auto const &[v1, v2] : itertools::zip(vi_out, vi_result)) {
    EXPECT_DOUBLE_EQ(v1, v2);
  }
  EXPECT_EQ(ui_from_vi(t_max, vi_out), std::vector<double>({6.0, 2.0, 1.0, 0.0}));
}

TEST(PlasmaUVWaroper, Mapping) { // NOLINT
  warper_plasma_uv_t w{10.0};
  double t_max = 10.0;
  std::vector<double> ui_times = {0.0, 2.0, 6.0, 1.0};

  std::vector<double> vi_out = w.li_from_ui(ui_times);
  std::vector<double> vi_result = {4.0, 4.0, 1.0, 1.0};

  EXPECT_EQ(vi_out, vi_result);

  std::vector<double> ui_sorted = {6.0, 2.0, 1.0, 0.0};
  EXPECT_EQ(w.ui_from_li(vi_out), ui_sorted);

  EXPECT_EQ(w.jacobian_reverse(vi_out), 1.0);
  EXPECT_EQ(w.jacobian_forward(ui_times), 1.0);
}

MAKE_MAIN; // NOLINT

// FIX OLD TESTS:

// double linear_function(double x) {
//   if (x > 1 || x < 0) {
//     return 0;
//   } else {
//     return x;
//   }
// }
// VI

// TEST(WarperPlasmaSimple, integrate_ansatz) { // NOLINT
//   double t_max       = 1.0;
//   int nr_grid_points = 3;
//   warper_product_1d_simple_interp_nearest_t warper(linear_function, t_max, nr_grid_points);

//   array<double, 1> f1_compare(nr_grid_points);
//   double delta = t_max / (nr_grid_points - 1);

//   triqs::clef::placeholder<0> i_;
//   f1_compare(i_) << delta * delta * i_ * i_;

//   for (auto const &[f1, f2] : itertools::zip(f1_compare, warper.f1_integrated.data())) { EXPECT_DOUBLE_EQ(f1, f2); }
// }