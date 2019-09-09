#include <keldy/impurity_oneband/compute_obs.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(ComputeObs, Initialize1) { // NOLINT
  model_param_t params;
  // TRIQS_PRINT(params.bath_type);
  compute_charge_Q computer(1, 10.0, params, 1000);
  // tests

  // run
  computer.run(1000);
  computer.run(1000);

  //
  std::cout << computer.reduce_result() << std::endl;
  std::cout << computer.reduce_nr_points_run() << std::endl;
  // get_warper

  // get_integrand
}

TEST(ComputeObs, Initialize2) { // NOLINT
  model_param_t params;
  // TRIQS_PRINT(params.bath_type);
  compute_charge_Q computer(4, 10.0, params, 1000);
  // tests

  // run
  computer.run(50000);
  //

  std::cout << computer.reduce_result() << std::endl;
  std::cout << computer.reduce_nr_points_run() << std::endl;
}

MAKE_MAIN; // NOLINT
