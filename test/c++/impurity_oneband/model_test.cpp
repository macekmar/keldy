#include <keldy/impurity_oneband/compute_obs.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(g0_model, Initialize) { // NOLINT
  model_param_t params;
  g0_model g0{params};
}

TEST(g0_keldysh_adaptor, Initialize) { // NOLINT
  model_param_t params;
  g0_model g0{params};
  g0_keldysh_contour_t g0_k{g0};
}

MAKE_MAIN; // NOLINT
