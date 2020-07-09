#!/usr/bin/env python

import unittest

from keldy import impurity_oneband as qqmc
from keldy import warpers
from keldy.visualization import integrand_warper_plot

import numpy as np
from mpi4py import MPI
import matplotlib
from matplotlib import pyplot as plt

class projection_usage(unittest.TestCase):

    def test_example(self):
        tmax = 20.
        params = {'beta': -1.,
                  'bias_V': 0.,
                  'eps_d': 0.,
                  'Gamma': 1.,
                  'alpha': 0.,
                  'half_bandwidth': 2.,
                  'bath_type': 'flatband',
                  'time_max': tmax,
                  'nr_time_points_gf': 10000,
                  'ft_method': 'analytic'}

        g0 = qqmc.G0Model(qqmc.G0ModelOmega(params), False)

        warper_train = warpers.WarperTrainT()
        warper_train.emplace_back(warpers.WarperPlasmaUvT(tmax))
        warper_train.emplace_back(warpers.make_product_1d_simple_exponential(15.0, 1.2))

        order = 3
        computer = qqmc.ComputeChargeQDirect(g0, tmax, order, 1e-15)
        computer.warper = warper_train

        warper_proj = warpers.WarperProjectionT(lambda u: np.abs(computer.evaluate_warped_integrand(u)[0]),
                                                order, 5000, 10000, 0.9, True)
        print(warper_proj.get_sigmas())

        warper_train.emplace_back(warper_proj)
        computer.warper = warper_train

        # ## plot the smoothing process
        # binner = warper_proj.get_xi(order - 1)
        # coord = binner.get_bin_coord()
        # plt.figure(0, (10, 8))
        # plt.plot(coord, binner.get_data(), '.', ms=1)
        # oplot(warper_proj.get_fi(order - 1))
        # plt.semilogy()
        # plt.show()


        for n in range(1, order + 1):
            integrand_warper_plot(computer, n, 0.1, tmax, nr_times=50)
            plt.show()

        computer.run(10000)

        print(computer.reduce_result())


if __name__ == '__main__':
    unittest.main()

