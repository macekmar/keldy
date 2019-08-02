#include <keldy/impurity_oneband/compute_obs.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(ComputeObs, Initialize1) { // NOLINT
  model_param_t params;
  compute_charge_Q computer(1, 5.0, params, 100);
  // tests

  // run
  computer.run(100);
  //
  dcomplex result = computer.reduce_result();
  std::cout << result << std::endl;
  // get_warper

  // get_integrand
}

TEST(ComputeObs, Initialize2) { // NOLINT
  model_param_t params;
  compute_charge_Q computer(2, 5.0, params, 100);
  // tests

  // run
  computer.run(100);
  //
  dcomplex result = computer.reduce_result();
  std::cout << result << std::endl;
  // get_warper

  // get_integrand
}

MAKE_MAIN; // NOLINT
