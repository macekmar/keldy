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


MAKE_MAIN; // NOLINT
