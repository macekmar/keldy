#include <keldy/impurity_oneband/compute_obs.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;


TEST(ComputeObs, Initialize2) { // NOLINT
  model_param_t params;
  compute_gf_kernel computer(params, 10.0, 4, "identity", 1000);

  computer.run(10);
  //

  // std::cout << computer.reduce_result() << std::endl;
  std::cout << computer.reduce_nr_points_run() << std::endl;



  
}


TEST(ComputeObs, Initialize3) { // NOLINT
  model_param_t params;
  compute_gf_kernel computer(params, 10.0, 4, "first_order", 1000);

  computer.run(10);

  auto result = computer.reduce_result();

  std::cout << result.get_bin_times() << std::endl;
  std::cout << result.get_bin_size() << std::endl;
  std::cout << result.get_nr_point_dropped() << std::endl;
  std::cout << result.get_nr_values()() << std::endl;
  std::cout << result.get_bin_times() << std::endl;

  //

  // std::cout << computer.reduce_result() << std::endl;
  std::cout << computer.reduce_nr_points_run() << std::endl;
}

MAKE_MAIN; // NOLINT
