"""
Utility functions for visualization of the integrand and/or warper.
"""

import numpy as np
from matplotlib import pyplot as plt

import keldy.impurity_oneband_module as keldy

def _gray_code(n):
    """
    Generates lists of n zeros and ones using the Gray code.
    """
    n = int(n)
    for i in range(0, 1<<n):
        gray = i^(i>>1)
        yield [int(c) for c in list("{0:0{1}b}".format(gray, n))]

def _ordered_axes(n):
    v_dir = np.zeros((n,), dtype=int)
    for i in range(n+1):
        v_dir[:i] = 1
        yield v_dir

def integrand_warper_plot(computer, order, d, t_max, nr_times=100, plot_all=False, axes=None, warper_prefactor=1.):
    """
    Plot the integrand weight and warper in absolute value along some axes in v space.

    A plot of |f(V + D)| is produced as a function of v, where:
      V = P({v, ..., v, 0, ..., 0})    with P a given permutation
      D = {d, d, ..., d}
    both are vectors of dimension `order`.
    f stands for the integrand or the warper from `computer`.

    `computer` is an integrator using the direct method
    `order` is the dimension of the integrand
    `d` is a positive real number
    `t_max` is the positive real t_max with which `computer` has been constructed
    `nr_times` is the number of points to plot
    `plot_all` is a boolean. If True all permutations P are displayed, if False only a subset of them (for clarity).
    `axes` (optionnal) a list of axes along which the integrand is plotted. Each axis is a list of 0s and 1s of length `order`. If provided, `plot_all` is ignored.
    `warper_prefactor` prefactor to warper in plot.
    """
    integrand = computer.get_integrand()
    warper = computer.get_warper()

    def list2str(l):
        s = ''
        for elt in l:
            if elt:
                s += 'v,'
            else:
                s += '0,'
        return s[:-1]

    def is_u_valid(u):
        return (u >= 0.).all() and (u <= t_max).all()

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    plt.plot([], [], '-k', label='integrand')
    plt.plot([], [], '--k', label='warper')

    x_arr = np.linspace(0., t_max, nr_times)

    if axes is None:
        if plot_all:
            v_dir_gen = _gray_code(order)
        else:
            v_dir_gen = _ordered_axes(order)
        v_dir_gen.next()
    else:
        v_dir_gen = axes

    for i, v_dir_int in enumerate(v_dir_gen):
        v_dir = np.array(v_dir_int, dtype=float)
        c = colors[i % len(colors)]

        v_arr = x_arr[:, None] * v_dir[None, :] + d * np.ones((1, order))
        u_arr = np.array([keldy.ui_from_vi(t_max, v) for v in v_arr])

        try: # for kernel
            values_v = [integrand(u)[0].sum_weights() if is_u_valid(u) else np.nan for u in u_arr]
        except AttributeError: # sum_weights not an attribute for non-kernel
            values_v = [np.abs(integrand(u)[0]) if is_u_valid(u) else np.nan for u in u_arr]

        plt.plot(x_arr, values_v, '-', c=c, label='V=({})'.format(list2str(v_dir_int)))

        values_v = [warper(u) if is_u_valid(u) else np.nan for u in u_arr]
        plt.plot(x_arr, warper_prefactor * np.abs(values_v), '--', c=c)

    plt.legend(loc=0)
    plt.xlabel(r'v')
    plt.ylabel(r'|f(V+D)|   (refer to docstring)')
    plt.title(r'Order n={}, shift d={}'.format(order, d))
    plt.semilogy()

    return x_arr

