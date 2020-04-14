#include "keldy/warpers/product_1d_simple.hpp"
#include <keldy/warpers/warpers.hpp>

#include <keldy/common.hpp>

#include <itertools/itertools.hpp>
#include <triqs/gfs.hpp>
#include <triqs/test_tools/arrays.hpp>
#include <triqs/test_tools/gfs.hpp>
#include <triqs/utility/macros.hpp>

using namespace keldy;
using namespace keldy::warpers;

using namespace triqs::arrays;

TEST(WarperTrain, ConstructionAddingWarpers) { // NOLINT
  warper_train_t w;

  double t_max = 10.0;
  w.warpers.emplace_back(warper_plasma_uv_t{t_max});

  // Copied from warper_plasma_uv_test
  std::vector<double> ui_times = {0.0, 2.0, 6.0, 1.0};

  std::vector<double> vi_out = w.li_from_ui(ui_times);
  std::vector<double> vi_result = {4.0, 4.0, 1.0, 1.0};

  EXPECT_EQ(vi_out, vi_result);

  std::vector<double> ui_sorted = {6.0, 2.0, 1.0, 0.0};
  EXPECT_EQ(w.ui_from_li(vi_out), ui_sorted);

  EXPECT_EQ(w.jacobian_reverse(vi_out), 1.0);
  EXPECT_EQ(w.jacobian_forward(ui_times), 1.0);

  // construct exponential warper
  double warper_scale = 1.0;
  w.warpers.emplace_back(warper_product_1d_simple_t{
     [warper_scale](double t) -> double { return std::exp(-(t / warper_scale)); },
     [warper_scale](double t) -> double { return warper_scale * (1 - std::exp(-t / warper_scale)); },
     [warper_scale](double l) -> double { return -warper_scale * std::log(1 - l / warper_scale); }, t_max, 100});

  // Make tests
}

MAKE_MAIN; // NOLINT
