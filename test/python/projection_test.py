#!/usr/bin/env python

import unittest

from keldy import impurity_oneband as qqmc
from keldy import warpers
from keldy.visualization import integrand_warper_plot, projection_warper_plot

import numpy as np
from mpi4py import MPI
import matplotlib
from matplotlib import pyplot as plt

class projection_usage(unittest.TestCase):

    def test_example(self):
        tmax = 20.
        order = 3
        params = {'beta': -1.,
                  'bias_V_left': 0.,
                  'bias_V_right': 0.,
                  'eps_d': 0.,
                  'Gamma': 1.,
                  'alpha': 0.,
                  'half_bandwidth': 2.,
                  'bath_type': 'flatband',
                  'time_max': order * tmax,
                  'nr_time_points_gf': 10000,
                  'ft_method': 'analytic'}

        g0 = qqmc.G0Model(qqmc.G0ModelOmega(params), False)

        warper_train = warpers.WarperTrainT()
        warper_train.emplace_back(warpers.WarperPlasmaUvT(tmax))
        warper_train.emplace_back(warpers.make_product_1d_simple_inverse(tmax, 1.2))

        computer = qqmc.ComputeChargeQDirect(g0, tmax, order, 1e-15)
        computer.warper = warper_train

        nr_bins = 500
        nr_samples = 5000
        warper_proj = warpers.WarperProjectionT(lambda l: computer.evaluate_warped_integrand(l, warper_train.size(), False),
                                                order, nr_bins, nr_samples, 0.1, True)
        print(warper_proj.get_sigmas())

        warper_train.emplace_back(warper_proj)
        computer.warper = warper_train

        ## plot the smoothing process
        projection_warper_plot(warper_proj, 1)

        for n in range(1, order + 1):
            integrand_warper_plot(computer, n, 0.1, tmax, nr_times=50)
            plt.show()

        computer.run(10000)

        print(computer.reduce_result())


if __name__ == '__main__':
    unittest.main()

