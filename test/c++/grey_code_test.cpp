#include <bitset>
#include <triqs/test_tools/gfs.hpp>

// using namespace keldy;
// using namespace keldy::anderson_re;

TEST(GreyCodeTest, Order3) { // NOLINT
  int order_n = 3;
  uint64_t nr_keldysh_configs = (uint64_t(1) << order_n);

  std::bitset<64> v = 0;

  for (int idx_kel = 0; idx_kel < nr_keldysh_configs; idx_kel++) {
    int nlc = (idx_kel < nr_keldysh_configs - 1 ? ffs(~idx_kel) : order_n) - 1;

    v[nlc] = !v[nlc];

    std::cout << "idx_kel: " << idx_kel << " ; nlc: " << nlc << " ; v: " << v << std::endl;
  }
  // model_param_t in_param;
  // model my_model(in_param);
  // std::cout << my_model.g0_lesser.mesh() << std::endl;
  // EXPECT_EQ(my_model.g0_lesser.mesh().size(), in_param.nr_time_points_gf);

  // EXPECT_EQ(my_model.g0_lesser.mesh().size(), in_param.nr_time_points_gf);
}

MAKE_MAIN // NOLINT
