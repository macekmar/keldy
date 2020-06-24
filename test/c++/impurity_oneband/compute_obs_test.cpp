#include "keldy/warpers/product_1d.hpp"
#include "keldy/warpers/product_1d_simple.hpp"
#include <keldy/impurity_oneband/compute_obs.hpp>
#include <triqs/test_tools/gfs.hpp>
#include <boost/math/constants/constants.hpp>

using namespace keldy;
using namespace keldy::impurity_oneband;

TEST(ComputeObs, Initialize1) { // NOLINT
  // model_param_t params;
  // // TRIQS_PRINT(params.bath_type);
  // compute_charge_Q_direct computer(params, 10.0, 1, 0.0, "first_order", 1000);
  // // tests

  // // run
  // computer.run(10);
  // computer.run(10);

  // //
  // std::cout << computer.reduce_result() << std::endl;
  // std::cout << computer.reduce_nr_points_run() << std::endl;
  // // get_warper

  // // get_integrand
}

TEST(ComputeObs, Initialize2) { // NOLINT
  // model_param_t params;
  // // TRIQS_PRINT(params.bath_type);
  // compute_charge_Q_direct computer(params, 10.0, 4, 0.0, "first_order", 1000);

  // // tests
  // // model_param_t params, double time, int order, std::string const &warper_function_name, int nr_sample_points_ansatz
  // // run
  // computer.run(10);
  // //

  // std::cout << computer.reduce_result() << std::endl;
  // std::cout << computer.reduce_nr_points_run() << std::endl;
}

TEST(ComputeObs, InitializeInverseSquare) { // NOLINT
  model_param_t params;
  double time = 10.0;
  compute_charge_Q_direct computer(params, time, 4, 0.0);

  computer.warper.emplace_back(warpers::warper_plasma_uv_t{time});
  computer.warper.emplace_back(warpers::make_product_1d_simple_inverse_square_nointerp(time, 1.5));

  computer.run(10);
}

TEST(ComputeObs, InitializeExponential) { // NOLINT
  model_param_t params;
  double time = 10.0;

  compute_charge_Q_direct computer(params, 10.0, 4, 0.0);

  computer.warper.emplace_back(warpers::warper_plasma_uv_t{time});
  computer.warper.emplace_back(warpers::make_product_1d_simple_exponential_nointerp(time, 1.5));

  computer.run(10);
}

TEST(ComputeObs, InitializeIdentity) { // NOLINT
  model_param_t params;
  double time = 10.0;

  compute_charge_Q_direct computer(params, time, 4, 0.0);
  computer.warper.emplace_back(warpers::warper_plasma_uv_t{time});

  computer.run(10);
}

TEST(ComputeObs, Initialize_current) { // NOLINT
  // model_param_t params;
  // // TRIQS_PRINT(params.bath_type);
  // compute_current_J_direct computer(params, 10.0, 4, 0.0, "first_order", 1000);

  // // tests
  // // model_param_t params, double time, int order, std::string const &warper_function_name, int nr_sample_points_ansatz
  // // run
  // computer.run(10);
  //

  // std::cout << computer.reduce_result() << std::endl;
  // std::cout << computer.reduce_nr_points_run() << std::endl;
}

TEST(ComputeChargeQDirect, FlatbandAnalyticCompare1) { // NOLINT

  auto i_to_the_n = std::vector<dcomplex>{1, 1_j, -1, -1_j};

  model_param_t params;
  params.beta = -1.0; // zero temperature
  params.bias_V = 0.0;
  params.eps_d = 0.0;
  params.Gamma = 1.0;
  params.time_max = 30.0;
  params.nr_time_points_gf = 10001;
  params.alpha = 0.0;
  params.bath_type = "flatband";
  params.ft_method = "analytic";

  // Series computed for epsilon_d=0, alpha=0, from Horvatic, Zlatic. J. Physique 46 (1985) 1459-1467.
  std::vector<double> HZ_charge_series = {0.5,
                                          -0.15915494309189533577,
                                          0.050660591821168885722,
                                          0.0046743460474992114797,
                                          -0.015882884762010443238,
                                          0.0071422698935926388204,
                                          0.0016484718800775438296};

  double tol = 1e-3;
  double t_max = 20.0;

  for (int order = 1; order <= 6; ++order) {
    std::cout << "Order " << order << std::endl;
    compute_charge_Q_direct computer(params, t_max, order, 0.0);

    computer.warper.emplace_back(warpers::warper_plasma_uv_t{t_max});
    computer.warper.emplace_back(warpers::make_product_1d_simple_inverse_square_nointerp(t_max, 1.0));

    computer.run(1e4);

    if (order >= 5) {
      tol = 1e-2;
    }

    EXPECT_COMPLEX_NEAR(-1_j * i_to_the_n[order % 4] * computer.reduce_result(), HZ_charge_series[order],
                        tol * std::abs(HZ_charge_series[order]));
  }
}

TEST(ComputeChargeQDirect, FlatbandAnalyticCompare2) { // NOLINT

  auto i_to_the_n = std::vector<dcomplex>{1, 1_j, -1, -1_j};

  model_param_t params;
  params.beta = -1.0; // zero temperature
  params.bias_V = 0.0;
  params.eps_d = 0.0;
  params.Gamma = 1.0;
  params.time_max = 30.0;
  params.nr_time_points_gf = 10001;
  params.alpha = 0.0;
  params.bath_type = "flatband";
  params.ft_method = "analytic";

  // Series computed for epsilon_d=0, alpha=0, from Horvatic, Zlatic. J. Physique 46 (1985) 1459-1467.
  std::vector<double> HZ_charge_series = {0.5,
                                          -0.15915494309189533577,
                                          0.050660591821168885722,
                                          0.0046743460474992114797,
                                          -0.015882884762010443238,
                                          0.0071422698935926388204,
                                          0.0016484718800775438296};

  double tol = 1e-2;
  double t_max = 20.0;

  for (int order = 1; order <= 6; ++order) {
    std::cout << "Order " << order << std::endl;
    compute_charge_Q_direct computer(params, t_max, order, 0.0);

    computer.warper.emplace_back(warpers::warper_plasma_uv_t{t_max});
    computer.warper.emplace_back(warpers::make_product_1d_simple_exponential_nointerp(t_max, 2.0));

    computer.run(5e4);

    //if (order >= 5) {
    //  tol = 1e-2;
    //}

    EXPECT_COMPLEX_NEAR(-1_j * i_to_the_n[order % 4] * computer.reduce_result(), HZ_charge_series[order],
                        tol * std::abs(HZ_charge_series[order]));
  }
}

TEST(ComputeObs, ValueCurrent1) { // NOLINT

  model_param_t params;
  params.beta = -1.0; // zero temperature
  params.bias_V = 0.6;
  params.eps_d = -0.5;
  params.Gamma = 1.0;
  params.time_max = 100.0;
  params.nr_time_points_gf = 5001;
  params.alpha = 0.;
  params.bath_type = "semicircle";
  params.ft_method = "contour";

  double t_max = 100.;
  compute_current_J_direct computer(params, t_max, 1, 0.0);

  computer.warper.emplace_back(warpers::warper_plasma_uv_t{t_max});
  computer.warper.emplace_back(warpers::make_product_1d_simple_inverse_square_nointerp(t_max, 0.5));

  computer.run(1e5);

  // exact result from analytical calculation
  EXPECT_NEAR(-2 * std::real(1_j * computer.reduce_result()), 0.041401946352562884, 1e-5);
}

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
  computer.warper.emplace_back(warpers::make_product_1d_inverse_cube_alternate(2, tmax, 1.0, 1e5));

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
