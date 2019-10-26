#include <keldy/common.hpp>
#include <keldy/impurity_oneband/wick_kernel.hpp>
#include <triqs/test_tools/gfs.hpp>
#include <triqs/arrays.hpp>

using triqs::arrays::array;
using triqs::arrays::range;

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(kernel_binner, Accumulation) { // NOLINT
  kernel_binner binner{0.0, 10.0, 4};

  // check bins are OK
  array<double, 1> times_expected(4);
  times_expected = {1.25, 3.75, 6.25, 8.75};
  EXPECT_ARRAY_EQ(times_expected, binner.get_bin_times());
  EXPECT_DOUBLE_EQ(2.5, binner.get_bin_size());

  binner.accumulate({1.5, up, forward},  1.0_j);
  binner.accumulate({6.3, up, backward}, 10.0_j);
  binner.accumulate({2.0, up, forward},  1.0);
  binner.accumulate({2.0, up, forward},  100.0_j);

  // boundary points
  binner.accumulate({0.0, up, backward}, 10.0);
  binner.accumulate({10.0, up, backward}, 100.0);

  // dropped values
  binner.accumulate({-1.0, up, forward}, 1000.0_j);
  binner.accumulate({11.0, up, forward}, 10'000.0_j);
  //TODO: check boundaries *between* bins

  // start testing
  EXPECT_EQ(2, binner.get_nr_point_dropped());

  array<dcomplex, 2> values_expected(4, 2);
  values_expected() = 0.0;
  values_expected(range(), 0) = array<dcomplex, 1>{1.0 + 101.0_j, 0.0_j, 0.0_j,  0.0_j        }; // forward
  values_expected(range(), 1) = array<dcomplex, 1>{10.0 + 0.0_j,  0.0_j, 10.0_j, 100.0 + 0.0_j}; // backward
  EXPECT_ARRAY_EQ(values_expected, binner.get_values());

  array<int, 2> nr_values_expected(4, 2);
  nr_values_expected() = 0;
  nr_values_expected(range(), 0) = array<int, 1>{3, 0, 0, 0};
  nr_values_expected(range(), 1) = array<int, 1>{1, 0, 1, 1};
  array<int, 2> nr_values(4, 2);
  nr_values() = binner.get_nr_values(); // cast to signed int
  EXPECT_ARRAY_EQ(nr_values_expected, nr_values);

  // product by scalar
  binner *= 2.0;
  EXPECT_ARRAY_EQ(2.0 * values_expected, binner.get_values());

}

TEST(kernel_binner, Accum_sparse_binner) { // NOLINT
  kernel_binner binner{0.0, 10.0, 4};
  sparse_kernel_binner sp_binner{{{{1.5, up, forward},  1.0_j},
                                  {{6.3, up, backward}, 10.0_j},
                                  {{2.0, up, forward},  1.0},
                                  {{2.0, up, forward},  100.0_j},
                                  {{0.0, up, backward}, 10.0},
                                  {{10.0, up, backward}, 100.0},
                                  {{-1.0, up, forward}, 1000.0_j},
                                  {{11.0, up, forward}, 10'000.0_j}}};

  binner += sp_binner;

  // start testing
  EXPECT_EQ(2, binner.get_nr_point_dropped());

  array<dcomplex, 2> values_expected(4, 2);
  values_expected() = 0.0;
  values_expected(range(), 0) = array<dcomplex, 1>{1.0 + 101.0_j, 0.0_j, 0.0_j,  0.0_j        }; // forward
  values_expected(range(), 1) = array<dcomplex, 1>{10.0 + 0.0_j,  0.0_j, 10.0_j, 100.0 + 0.0_j}; // backward
  EXPECT_ARRAY_EQ(values_expected, binner.get_values());

  array<int, 2> nr_values_expected(4, 2);
  nr_values_expected() = 0;
  nr_values_expected(range(), 0) = array<int, 1>{3, 0, 0, 0};
  nr_values_expected(range(), 1) = array<int, 1>{1, 0, 1, 1};
  array<int, 2> nr_values(4, 2);
  nr_values() = binner.get_nr_values(); // cast to signed int
  EXPECT_ARRAY_EQ(nr_values_expected, nr_values);
}

TEST(sparse_kernel_binner, weight) { // NOLINT
  sparse_kernel_binner sp_binner{{{{1.5, up, forward},  1.0_j},
                                 {{6.3, up, backward}, 10.0_j},
                                 {{2.0, up, forward},  1.0},
                                 {{2.0, up, forward},  100.0_j}}};

  EXPECT_DOUBLE_EQ(112.0, sp_binner.sum_weights());
}

TEST(sparse_kernel_binner, product_scalar) { // NOLINT
  sparse_kernel_binner sp_binner{{{{1.5, up, forward},  1.0_j},
                                 {{6.3, up, backward}, 10.0_j},
                                 {{2.0, up, forward},  1.0},
                                 {{2.0, up, forward},  100.0_j}}};
  sp_binner *= 2.0;

  EXPECT_DOUBLE_EQ(2.0 * 112.0, sp_binner.sum_weights());
  // TODO: check values
}

MAKE_MAIN; // NOLINT
