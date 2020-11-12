#include "keldy/common.hpp"
#include <keldy/impurity_oneband/compute_obs.hpp>
#include <triqs/arrays.hpp>

#include <triqs/test_tools/gfs.hpp>
#include <algorithm>

using namespace keldy;
using namespace keldy::impurity_oneband;
namespace nda = triqs::arrays;

TEST(ChiFunction, Initialization) { // NOLINT
  model_param_t const params;
  g0_model const g0{g0_model_omega{params}, false};
  g0_keldysh_contour_t const g0_k{g0};

  nda::array<double, 1> omegas = {0.0, 0.5, 1.0, 1.5, 2.0};
  chi_function_t const chi(g0.model_omega, omegas, 10.0, backward);

  gf_index_t A(2.0, up, forward);
  gf_index_t B(5.0, down, backward);

  std::cout << chi(A, B) << std::endl;
}

TEST(SpinPlusSpinMinus, Initialization) { //NOLINT
  model_param_t const params;
  nda::array<double, 1> omegas = {0.0, 0.5, 1.0, 1.5, 2.0};

  compute_spin_plus_spin_minus_freq computer(params, 10.0, omegas, 3);
  auto integrand = computer.get_integrand();

  std::cout << integrand({2.0, 5.0, 7.0}) << std::endl;
  std::cout << computer.evaluate_warped_integrand({2.0, 5.0, 7.0}) << std::endl;

  computer.run(10);

  std::cout << computer.reduce_result() << std::endl;
  std::cout << computer.reduce_nr_points_run() << std::endl;
  std::cout << computer.reduce_nr_points_in_domain() << std::endl;
}

MAKE_MAIN // NOLINT
