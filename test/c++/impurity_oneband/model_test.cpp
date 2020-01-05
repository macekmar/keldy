#include "keldy/common.hpp"
#include "keldy/impurity_oneband/model.hpp"
#include "keldy/contour_integral.hpp"
#include <gtest/gtest.h>
#include <itertools/itertools.hpp>
#include <keldy/impurity_oneband/compute_obs.hpp>
//#include <keldy/impurity_oneband/model.hpp>
#include <triqs/test_tools/gfs.hpp>
//#include <cmath>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/constants/constants.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(g0_keldysh_adaptor, Initialize) { // NOLINT
  model_param_t params;
  g0_model g0{params, true};
  g0_keldysh_contour_t g0_k{g0};
}

MAKE_MAIN; // NOLINT
