#include <gtest/gtest.h>
#include <keldy/qrng.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy;

TEST(RNGConstructors, SobolUnshifted) { // NOLINT
  int dimension = 3;
  int rng_state_seed = 0;
  auto rng = sobol(dimension, rng_state_seed);

  auto output = rng();

  EXPECT_EQ(output.size(), dimension);
  for (auto o : output) {
    std::cout << o << ",";
  }
  std::cout << std::endl;
}

TEST(RNGConstructors, MT) { // NOLINT
  int dimension = 3;
  int rng_state_seed = 0;
  auto rng = boost_rng_wrapper_t<boost::random::mt19937_64>(dimension, rng_state_seed);

  auto output = rng();

  EXPECT_EQ(output.size(), dimension);
  for (auto o : output) {
    std::cout << o << ",";
  }
  std::cout << std::endl;
}

TEST(RNGConstructors, Randlux) { // NOLINT
  int dimension = 3;
  int rng_state_seed = 0;
  auto rng = boost_rng_wrapper_t<boost::random::ranlux48>(dimension, rng_state_seed);

  auto output = rng();

  EXPECT_EQ(output.size(), dimension);
  for (auto o : output) {
    std::cout << o << ",";
  }
  std::cout << std::endl;
}

TEST(RNGConstructors, Taus) { // NOLINT
  int dimension = 3;
  int rng_state_seed = 0;
  auto rng = boost_rng_wrapper_t<boost::random::taus88>(dimension, rng_state_seed);

  auto output = rng();

  EXPECT_EQ(output.size(), dimension);
  for (auto o : output) {
    std::cout << o << ",";
  }
  std::cout << std::endl;
}

MAKE_MAIN // NOLINT
