#include "keldy/binner.hpp"
#include "tests_std_containers.hpp"
#include <triqs/test_tools/gfs.hpp>
#include <gtest/gtest.h>

using namespace keldy::binner;
using namespace triqs::arrays;

TEST(Binner, Construction) { // NOLINT
  auto bin10 = binner_t<1>({std::make_tuple(-2.5, 5.0, 10)});
  bin10.accumulate(1_j, 0.5);
  bin10.accumulate(1_j, 0.5);

  auto bin11 = binner_t<1, 1>({std::make_tuple(-2.5, 5.0, 10)}, {3});
  bin11.accumulate(1_j, 0.5, 3);
  bin11.accumulate(1_j, 0.5, 1);

  auto bin12 = binner_t<1, 2>({std::make_tuple(-2.5, 5.0, 10)}, {3, 6});
  bin12.accumulate(1_j, 0.5, 3, 11);
  bin12.accumulate(1_j, 0.5, 1, 2);

  auto bin20 = binner_t<2>({std::make_tuple(-2.5, 5.0, 10), std::make_tuple(0., 10., 5)});
  bin20.accumulate(1_j, 0.5, 3.);
  bin20.accumulate(1_j, 0.5, 13.);

  auto bin21 = binner_t<2, 1>({std::make_tuple(-2.5, 5.0, 10), std::make_tuple(0., 10., 5)}, {3});
  bin21.accumulate(1_j, 0.5, 3., 0);
  bin21.accumulate(1_j, 0.5, 13., 5);
}

TEST(Binner, accumulate) { // NOLINT
  binner_t<1, 1> binner({std::make_tuple(0.0, 10.0, 4)}, {2});

  // check bins are OK
  array<double, 1> times_expected(4);
  times_expected = {1.25, 3.75, 6.25, 8.75};
  EXPECT_ARRAY_EQ(times_expected, binner.get_bin_coord());
  EXPECT_DOUBLE_EQ(2.5, binner.get_bin_size());

  binner.accumulate(1.0_j, 1.5, 0);
  binner.accumulate(10.0_j, 6.3, 1);
  binner.accumulate(1.0, 2.0, 0);
  binner.accumulate(100.0_j, 2.0, 0);

  // boundary points
  binner.accumulate(10.0, 0.0, 1);
  binner.accumulate(100.0, 10.0, 1);
  binner.accumulate(1000.0, 2.5, 1);
  binner.accumulate(10'000.0, 5.0, 1);
  binner.accumulate(100'000.0, 7.5, 1);

  // dropped values
  binner.accumulate(1000.0_j, -1.0, 0);
  binner.accumulate(10'000.0_j, 11.0, 0);

  // start testing
  EXPECT_EQ(2, binner.get_nr_values_dropped());

  array<dcomplex, 2> values_expected(4, 2);
  values_expected() = 0.0;
  values_expected(range(), 0) = array<dcomplex, 1>{1.0 + 101.0_j, 0.0_j, 0.0_j, 0.0_j}; // forward
  values_expected(range(), 1) =
     array<dcomplex, 1>{10.0 + 0.0_j, 1000.0 + 0.0_j, 10'000.0 + 10.0_j, 100'100.0 + 0.0_j}; // backward
  EXPECT_ARRAY_EQ(values_expected, binner.get_data());

  array<long, 2> nr_values_expected(4, 2);
  nr_values_expected() = 0;
  nr_values_expected(range(), 0) = array<int, 1>{3, 0, 0, 0};
  nr_values_expected(range(), 1) = array<int, 1>{1, 1, 2, 2};
  array<long, 2> nr_values = binner.get_nr_values_added(); // cast to signed long so that difference can be taken
  EXPECT_ARRAY_EQ(nr_values_expected, nr_values);

  // product by scalar
  binner *= 2.0;
  EXPECT_ARRAY_EQ(2.0 * values_expected, binner.get_data());
}

TEST(SparseBinner, Construction) { // NOLINT
  auto bin10 = sparse_binner_t<1>();
  bin10.accumulate(1_j, 0.5);

  auto bin11 = sparse_binner_t<1, 1>();
  bin11.accumulate(1_j, 0.5, 3);

  auto bin12 = sparse_binner_t<1, 2>();
  bin12.accumulate(1_j, 0.5, 3, 11);

  auto bin20 = sparse_binner_t<2>();
  bin20.accumulate(1_j, 0.5, 3.);

  auto bin21 = sparse_binner_t<2, 1>();
  bin21.accumulate(1_j, 0.5, 3., 0);
}

TEST(SparseBinner, accumulate) { //NOLINT
  sparse_binner_t<1, 2> sp_binner;
  sp_binner.accumulate(1_j, 1.5, 0, 3);
  sp_binner.accumulate(10_j, 1.5, 2, 3);
  sp_binner.accumulate(100_j, 1.5, 2, 3);
  sp_binner.accumulate(1000_j, 1.50001, 2, 3);

  auto &data = sp_binner.data;

  EXPECT_EQ(data.size(), 3);
  EXPECT_EQ(data[0].second, 1_j);
  EXPECT_EQ(data[1].second, 110_j);
  EXPECT_EQ(data[2].second, 1000_j);

  EXPECT_TRUE(are_iterable_eq(data[0].first.first, {1.5}));
  EXPECT_TRUE(are_iterable_eq(data[0].first.second, {0, 3}));
  EXPECT_TRUE(are_iterable_eq(data[1].first.first, {1.5}));
  EXPECT_TRUE(are_iterable_eq(data[1].first.second, {2, 3}));
  EXPECT_TRUE(are_iterable_eq(data[2].first.first, {1.50001}));
  EXPECT_TRUE(are_iterable_eq(data[2].first.second, {2, 3}));
}

TEST(SparseBinner, weight) { // NOLINT
  sparse_binner_t<1, 1> sp_binner;
  sp_binner.accumulate(1.0_j, 1.5, 0);
  sp_binner.accumulate(10.0_j, 6.3, 1);
  sp_binner.accumulate(100.0_j, 2.0, 0);
  sp_binner.accumulate(1000.0_j, 2.0, 0);

  EXPECT_DOUBLE_EQ(1111.0, sp_binner.sum_moduli());
}

TEST(SparseBinner, product_scalar) { // NOLINT
  sparse_binner_t<1, 1> sp_binner;
  sp_binner.accumulate(1.0_j, 1.5, 0);
  sp_binner.accumulate(10.0_j, 6.3, 1);
  sp_binner.accumulate(100.0_j, 2.0, 0);
  sp_binner.accumulate(1000.0_j, 2.0, 0);

  sp_binner *= 2.0;

  EXPECT_DOUBLE_EQ(2.0 * 1111.0, sp_binner.sum_moduli());
}

TEST(BinnerAndSparseBinner, accumulate) { // NOLINT
  binner_t<1, 1> binner({std::make_tuple(0.0, 10.0, 4)}, {2});
  sparse_binner_t<1, 1> sp_binner;
  sp_binner.accumulate(1.0_j, 1.5, 0);
  sp_binner.accumulate(10.0_j, 6.3, 1);
  sp_binner.accumulate(1.0, 2.0, 0);     // same coord
  sp_binner.accumulate(100.0_j, 2.0, 0); // same coord
  sp_binner.accumulate(10.0, 0.0, 1);
  sp_binner.accumulate(100.0, 10.0, 1);
  sp_binner.accumulate(1000.0, 2.5, 1);
  sp_binner.accumulate(10'000.0, 5.0, 1);
  sp_binner.accumulate(100'000.0, 7.5, 1);
  sp_binner.accumulate(1000.0_j, -1.0, 0);
  sp_binner.accumulate(10'000.0_j, 11.0, 0);

  binner += sp_binner;

  // start testing
  EXPECT_EQ(2, binner.get_nr_values_dropped());

  array<dcomplex, 2> values_expected(4, 2);
  values_expected() = 0.0;
  values_expected(range(), 0) = array<dcomplex, 1>{1.0 + 101.0_j, 0.0_j, 0.0_j, 0.0_j}; // forward
  values_expected(range(), 1) =
     array<dcomplex, 1>{10.0 + 0.0_j, 1000.0 + 0.0_j, 10'000.0 + 10.0_j, 100'100.0 + 0.0_j}; // backward
  EXPECT_ARRAY_EQ(values_expected, binner.get_data());

  array<long, 2> nr_values_expected(4, 2);
  nr_values_expected() = 0;
  nr_values_expected(range(), 0) = array<int, 1>{2, 0, 0, 0};
  nr_values_expected(range(), 1) = array<int, 1>{1, 1, 2, 2};
  array<long, 2> nr_values = binner.get_nr_values_added(); // cast to signed long so that difference can be taken
  EXPECT_ARRAY_EQ(nr_values_expected, nr_values);
}

TEST(BinnerAndSparseBinner, accumulate_twice) { // NOLINT
  binner_t<1, 1> binner({std::make_tuple(0.0, 10.0, 4)}, {2});
  sparse_binner_t<1, 1> sp_binner_A;
  sp_binner_A.accumulate(1.0_j, 1.5, 0);
  sp_binner_A.accumulate(10.0_j, 6.3, 1);
  sp_binner_A.accumulate(1.0, 2.0, 0);     // same coord
  sp_binner_A.accumulate(100.0_j, 2.0, 0); // same coord
  sp_binner_A.accumulate(10.0, 0.0, 1);
  sp_binner_A.accumulate(100.0, 10.0, 1);

  binner += sp_binner_A;

  sparse_binner_t<1, 1> sp_binner_B;
  sp_binner_B.accumulate(1000.0, 2.5, 1);
  sp_binner_B.accumulate(10'000.0, 5.0, 1);
  sp_binner_B.accumulate(100'000.0, 7.5, 1);
  sp_binner_B.accumulate(1000.0_j, -1.0, 0);
  sp_binner_B.accumulate(10'000.0_j, 11.0, 0);

  binner += sp_binner_B;

  // start testing
  EXPECT_EQ(2, binner.get_nr_values_dropped());

  array<dcomplex, 2> values_expected(4, 2);
  values_expected() = 0.0;
  values_expected(range(), 0) = array<dcomplex, 1>{1.0 + 101.0_j, 0.0_j, 0.0_j, 0.0_j}; // forward
  values_expected(range(), 1) =
     array<dcomplex, 1>{10.0 + 0.0_j, 1000.0 + 0.0_j, 10'000.0 + 10.0_j, 100'100.0 + 0.0_j}; // backward
  EXPECT_ARRAY_EQ(values_expected, binner.get_data());

  array<long, 2> nr_values_expected(4, 2);
  nr_values_expected() = 0;
  nr_values_expected(range(), 0) = array<int, 1>{2, 0, 0, 0};
  nr_values_expected(range(), 1) = array<int, 1>{1, 1, 2, 2};
  array<long, 2> nr_values = binner.get_nr_values_added(); // cast to signed long so that difference can be taken
  EXPECT_ARRAY_EQ(nr_values_expected, nr_values);
}

MAKE_MAIN; // NOLINT
