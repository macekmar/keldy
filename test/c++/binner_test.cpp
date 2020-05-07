#include "keldy/binner.hpp"
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
  bin10.append(1_j, 0.5);

  auto bin11 = sparse_binner_t<1, 1>();
  bin11.accumulate(1_j, 0.5, 3);
  bin11.append(1_j, 0.5, 1);

  auto bin12 = sparse_binner_t<1, 2>();
  bin12.accumulate(1_j, 0.5, 3, 11);
  bin12.append(1_j, 0.5, 1, 2);

  auto bin20 = sparse_binner_t<2>();
  bin20.accumulate(1_j, 0.5, 3.);
  bin20.append(1_j, 0.5, 13.);

  auto bin21 = sparse_binner_t<2, 1>();
  bin21.accumulate(1_j, 0.5, 3., 0);
  bin21.append(1_j, 0.5, 13., 5);
}

TEST(SparseBinner, accumulate) { //NOLINT
  sparse_binner_t<1, 2> sp_binner;
  sp_binner.accumulate(1_j, 1.5, 0, 3);
  sp_binner.accumulate(10_j, 1.5, 2, 3);
  sp_binner.accumulate(100_j, 1.5, 2, 3);
  sp_binner.accumulate(1000_j, 1.50001, 2, 3);

  EXPECT_EQ(sp_binner.data.size(), 3);
  EXPECT_EQ(sp_binner.data[0].second, 1_j);
  EXPECT_EQ(sp_binner.data[1].second, 110_j);
  EXPECT_EQ(sp_binner.data[2].second, 1000_j);

  //TODO: check coordinates
}

TEST(SparseBinner, append) { //NOLINT
  sparse_binner_t<1, 2> sp_binner;
  sp_binner.append(1_j, 1.5, 0, 3);
  sp_binner.append(10_j, 1.5, 2, 3);
  sp_binner.append(100_j, 1.5, 2, 3);
  sp_binner.append(1000_j, 1.50001, 2, 3);

  EXPECT_EQ(sp_binner.data.size(), 4);
  EXPECT_EQ(sp_binner.data[0].second, 1_j);
  EXPECT_EQ(sp_binner.data[1].second, 10_j);
  EXPECT_EQ(sp_binner.data[2].second, 100_j);
  EXPECT_EQ(sp_binner.data[3].second, 1000_j);

  //TODO: check coordinates
}

TEST(SparseBinner, weight) { // NOLINT
  sparse_binner_t<1, 1> sp_binner;
  sp_binner.append(1.0_j, 1.5, 0);
  sp_binner.append(10.0_j, 6.3, 1);
  sp_binner.append(1.0, 2.0, 0);
  sp_binner.append(100.0_j, 2.0, 0);

  //std::cout << sp_binner.data.size() << std::endl;
  //std::cout << sp_binner.data[2].first << std::endl;
  //std::cout << sp_binner.data[2].second << std::endl;

  EXPECT_DOUBLE_EQ(112.0, sp_binner.sum_moduli());
}

TEST(SparseBinner, product_scalar) { // NOLINT
  sparse_binner_t<1, 1> sp_binner;
  sp_binner.append(1.0_j, 1.5, 0);
  sp_binner.append(10.0_j, 6.3, 1);
  sp_binner.append(1.0, 2.0, 0);
  sp_binner.append(100.0_j, 2.0, 0);

  sp_binner *= 2.0;

  EXPECT_DOUBLE_EQ(2.0 * 112.0, sp_binner.sum_moduli());
  EXPECT_EQ(sp_binner.data[0].second, 2.0_j);
  EXPECT_EQ(sp_binner.data[1].second, 20.0_j);
  EXPECT_EQ(sp_binner.data[2].second, 2.0);
  EXPECT_EQ(sp_binner.data[3].second, 200.0_j);
}

TEST(SparseBinner, equality) { //NOLINT
  sparse_binner_t<2, 1> sp_bin_A;
  sp_bin_A.append(1.0, 3.5, -0.7, 8);
  sp_bin_A.append(10.0, 3.5, -0.7, 8); // same coord
  sp_bin_A.append(100.0, 1.0, -0.7, 5);
  sp_bin_A.append(1000.0, 2.0, -1.0, 4);

  sparse_binner_t<2, 1> sp_bin_B; // same but populated in different order
  sp_bin_B.append(1000.0, 2.0, -1.0, 4);
  sp_bin_B.append(1.0, 3.5, -0.7, 8);
  sp_bin_B.append(100.0, 1.0, -0.7, 5);
  sp_bin_B.append(10.0, 3.5, -0.7, 8); // same coord

  EXPECT_TRUE(sp_bin_A == sp_bin_B);

  sp_bin_A.append(1.0_j, 4.0, 1.0, 2);
}

TEST(SparseBinner, non_equality_size) { //NOLINT
  sparse_binner_t<2, 1> sp_bin_A;
  sp_bin_A.append(1.0, 3.5, -0.7, 8);
  sp_bin_A.append(100.0, 1.0, -0.7, 5);
  sp_bin_A.append(1.0_j, 4.5, 1.0, 2);

  sparse_binner_t<2, 1> sp_bin_B;
  sp_bin_B.append(1.0, 3.5, -0.7, 8);
  sp_bin_B.append(100.0, 1.0, -0.7, 5);

  EXPECT_FALSE(sp_bin_A == sp_bin_B);
}

TEST(SparseBinner, non_equality_coord) { //NOLINT
  sparse_binner_t<2, 1> sp_bin_A;
  sp_bin_A.append(1.0, 3.5, -0.7, 8);
  sp_bin_A.append(100.0, 1.0, -0.7, 5);

  sp_bin_A.append(1.0_j, 4.5, 1.0, 2);

  sparse_binner_t<2, 1> sp_bin_B;
  sp_bin_B.append(1.0, 3.5, -0.7, 8);
  sp_bin_B.append(100.0, 1.0, -0.7, 5);

  sp_bin_B.append(1.0_j, 4.0, 1.0, 2);

  EXPECT_FALSE(sp_bin_A == sp_bin_B);
}

TEST(SparseBinner, non_equality_value) { //NOLINT
  sparse_binner_t<2, 1> sp_bin_A;
  sp_bin_A.append(1.0, 3.5, -0.7, 8);
  sp_bin_A.append(100.0, 1.0, -0.7, 5);

  sp_bin_A.append(1.0_j, 4.0, 1.0, 2);

  sparse_binner_t<2, 1> sp_bin_B;
  sp_bin_B.append(1.0, 3.5, -0.7, 8);
  sp_bin_B.append(100.0, 1.0, -0.7, 5);

  sp_bin_B.append(10.0_j, 4.0, 1.0, 2);

  EXPECT_FALSE(sp_bin_A == sp_bin_B);
}

/// TODO: test sparse_binner_t sort, test if it is stable

TEST(BinnerAndSparseBinner, accumulate) { // NOLINT
  binner_t<1, 1> binner({std::make_tuple(0.0, 10.0, 4)}, {2});
  sparse_binner_t<1, 1> sp_binner;
  sp_binner.append(1.0_j, 1.5, 0);
  sp_binner.append(10.0_j, 6.3, 1);
  sp_binner.append(1.0, 2.0, 0);
  sp_binner.append(100.0_j, 2.0, 0);
  sp_binner.append(10.0, 0.0, 1);
  sp_binner.append(100.0, 10.0, 1);
  sp_binner.append(1000.0, 2.5, 1);
  sp_binner.append(10'000.0, 5.0, 1);
  sp_binner.append(100'000.0, 7.5, 1);
  sp_binner.append(1000.0_j, -1.0, 0);
  sp_binner.append(10'000.0_j, 11.0, 0);

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
  nr_values_expected(range(), 0) = array<int, 1>{3, 0, 0, 0};
  nr_values_expected(range(), 1) = array<int, 1>{1, 1, 2, 2};
  array<long, 2> nr_values = binner.get_nr_values_added(); // cast to signed long so that difference can be taken
  EXPECT_ARRAY_EQ(nr_values_expected, nr_values);
}

TEST(BinnerAndSparseBinner, accumulate_twice) { // NOLINT
  binner_t<1, 1> binner({std::make_tuple(0.0, 10.0, 4)}, {2});
  sparse_binner_t<1, 1> sp_binner_A;
  sp_binner_A.append(1.0_j, 1.5, 0);
  sp_binner_A.append(10.0_j, 6.3, 1);
  sp_binner_A.append(1.0, 2.0, 0);
  sp_binner_A.append(100.0_j, 2.0, 0);
  sp_binner_A.append(10.0, 0.0, 1);
  sp_binner_A.append(100.0, 10.0, 1);

  binner += sp_binner_A;

  sparse_binner_t<1, 1> sp_binner_B;
  sp_binner_B.append(1000.0, 2.5, 1);
  sp_binner_B.append(10'000.0, 5.0, 1);
  sp_binner_B.append(100'000.0, 7.5, 1);
  sp_binner_B.append(1000.0_j, -1.0, 0);
  sp_binner_B.append(10'000.0_j, 11.0, 0);

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
  nr_values_expected(range(), 0) = array<int, 1>{3, 0, 0, 0};
  nr_values_expected(range(), 1) = array<int, 1>{1, 1, 2, 2};
  array<long, 2> nr_values = binner.get_nr_values_added(); // cast to signed long so that difference can be taken
  EXPECT_ARRAY_EQ(nr_values_expected, nr_values);
}

MAKE_MAIN; // NOLINT
