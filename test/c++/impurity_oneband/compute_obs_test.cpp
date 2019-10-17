#include <keldy/impurity_oneband/compute_obs.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(ComputeObs, Initialize1) { // NOLINT
  model_param_t params;
  // TRIQS_PRINT(params.bath_type);
  compute_charge_Q_direct computer(params, 10.0, 1, "first_order", 1000);
  // tests

  // run
  computer.run(10);
  computer.run(10);

  //
  std::cout << computer.reduce_result() << std::endl;
  std::cout << computer.reduce_nr_points_run() << std::endl;
  // get_warper

  // get_integrand
}

TEST(ComputeObs, Initialize2) { // NOLINT
  model_param_t params;
  // TRIQS_PRINT(params.bath_type);
  compute_charge_Q_direct computer(params, 10.0, 4, "first_order", 1000);
  
  // tests
// model_param_t params, double time, int order, std::string const &warper_function_name, int nr_sample_points_ansatz
  // run
  computer.run(10);
  //

  std::cout << computer.reduce_result() << std::endl;
  std::cout << computer.reduce_nr_points_run() << std::endl;
}

MAKE_MAIN; // NOLINT
