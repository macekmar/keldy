#include "keldy/warpers/product_1d.hpp"
#include "keldy/warpers/product_1d_simple.hpp"
#include <keldy/impurity_oneband/compute_obs.hpp>
#include <triqs/test_tools/gfs.hpp>
#include <boost/math/constants/constants.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(ComputeObs, ValueKernel) { // NOLINT
  double const tmax = 100.;

  model_param_t params;
  params.beta = -1.0; // zero temperature
  params.bias_V = 0.;
  params.eps_d = 0.;
  params.Gamma = 0.5;
  params.alpha = 0.5;
  params.half_bandwidth = 6.;
  params.time_max = tmax;
  params.nr_time_points_gf = 5001;
  params.bath_type = "semicircle";
  params.ft_method = "contour";

  auto const g0 = g0_model(g0_model_omega(params), false);

  int const nr_bins = 10;
  compute_gf_kernel computer(g0, tmax, 2, nr_bins);

  computer.warper.emplace_back(warpers::warper_plasma_uv_t{tmax});
  computer.warper.emplace_back(warpers::make_product_1d_inverse_cube_alternate_interp_hybrid(2, tmax, 1.0, 1e5));

  computer.run(1e6);
  auto const binner = computer.reduce_result();
  auto const result = make_array(binner.get_data() / binner.get_bin_size());

  /// lesser is the conjugate of greater (particle-hole symmetry)
  EXPECT_LT(max_element(abs(result(range(), 0) - conj(result(range(), 1)))), 1e-12);

  /// previous calculation (job #524060, 26 March 2020, probably on commit 10319ee) gives the following, with accuracy 2e-8.
  auto const expected_result =
     array<dcomplex, 2>{{-3.67742471e-09 + 1.63164724e-04_j, -3.67742341e-09 - 1.63164724e-04_j},
                        {5.50958865e-09 + 2.09202817e-05_j, 5.50959042e-09 - 2.09202817e-05_j},
                        {-7.41242821e-09 + 2.68650657e-05_j, -7.41242621e-09 - 2.68650657e-05_j},
                        {9.25484502e-09 + 3.59362761e-05_j, 9.25484734e-09 - 3.59362761e-05_j},
                        {-1.08196801e-08 + 5.05717360e-05_j, -1.08196774e-08 - 5.05717360e-05_j},
                        {1.16796193e-08 + 7.66168036e-05_j, 1.16796227e-08 - 7.66168036e-05_j},
                        {-1.50691787e-08 + 1.30292898e-04_j, -1.50691743e-08 - 1.30292898e-04_j},
                        {-9.34236787e-07 + 2.77177129e-04_j, -9.34236781e-07 - 2.77177129e-04_j},
                        {-2.00377184e-04 + 1.06584119e-03_j, -2.00377184e-04 - 1.06584119e-03_j},
                        {2.01317477e-04 - 1.85006394e-03_j, 2.01317477e-04 + 1.85006394e-03_j}};

  std::cout << max_element(abs(result - expected_result)) << std::endl;
  EXPECT_LT(max_element(abs(result - expected_result)), 1e-6);

  std::cout << result << std::endl;
}

MAKE_MAIN; // NOLINT
