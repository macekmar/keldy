#include <keldy/common.hpp>
#include <keldy/impurity_oneband/wick_direct.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;


TEST(g0_keldysh_adaptor, Initialize) { // NOLINT
  model_param_t params;
  g0_model g0{params};
  g0_keldysh_contour_t g0_k{g0};
  integrand_g_direct integrand{g0_k, gf_index_t{5.0, up, forward}, gf_index_t{5.0, up, forward}};

  integrand(std::vector<double>{1.0});

  integrand(std::vector<double>{1.0, 2.0});

}

MAKE_MAIN; // NOLINT
