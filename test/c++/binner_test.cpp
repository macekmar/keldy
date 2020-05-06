#include "keldy/binner.hpp"
#include <triqs/test_tools/gfs.hpp>
#include <gtest/gtest.h>

using namespace keldy::binner;
using namespace triqs::arrays;

TEST(Binner, Construction_1) { // NOLINT

  auto bin = binner_t<1, 2>({std::make_tuple(-2.5, 5.0, 10)}, {4, 2});

  bin.accumulate(1_j, 0.5, 3, 0);
  bin.accumulate(1_j, 0.7, 2, 1);
  bin.accumulate(1_j, -6.0, 2, 1); // out of binner
  bin.accumulate(1_j, 3.0, 8, 1);  // out of binner

  EXPECT_EQ(bin.get_nr_values_dropped(), 2);

  array<dcomplex, 3> data = bin.get_data();
  EXPECT_EQ(first_dim(data), 10);
  EXPECT_EQ(second_dim(data), 4);
  EXPECT_EQ(third_dim(data), 2);
  std::cout << fourth_dim(data) << std::endl; // returns a number ?????
  std::cout << get_shape(data) << std::endl;
  std::cout << data << std::endl;

  data(0, 0, 0) = 10_j;
  EXPECT_NE(data(0, 0, 0), bin.get_data()(0, 0, 0));

  array<unsigned long, 3> nr_values = bin.get_nr_values_added();
  EXPECT_EQ(first_dim(nr_values), 10);
  EXPECT_EQ(second_dim(nr_values), 4);
  EXPECT_EQ(third_dim(nr_values), 2);
  std::cout << nr_values << std::endl;
}

TEST(Binner, Construction_2) { // NOLINT
  auto bin = binner_t<2>({std::make_tuple(-2.5, 5.0, 10), std::make_tuple(0., 10., 5)});

  bin.accumulate(1_j, 0.5, 3.);
  bin.accumulate(1_j, 0.5, 13.);
}

TEST(Binner, sparse_binner) { // NOLINT
  sparse_binner_t<1, 2> sp_bin;

  sp_bin.accumulate(1_j, 0.5, 3, 0);
  sp_bin.accumulate(1_j, 0.7, 2, 1);
  sp_bin.accumulate(1_j, -6.0, 2, 1); // out of binner
  sp_bin.accumulate(1_j, 3.0, 8, 1);  // out of binner

  auto res = sp_bin.operator*=(5.);

  auto bin = binner_t<1, 2>({std::make_tuple(-2.5, 5.0, 10)}, {4, 2});

  bin += sp_bin;

  std::cout << bin.get_data() << std::endl;
}

MAKE_MAIN; // NOLINT
