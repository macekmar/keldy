"""
Utility functions to extract retarded kernel and Green function.
"""
import numpy as _np
from numpy import fft as _fft


# stolen from scipy.signal.fftconvolve
def _next_regular(target):
    """
    Find the next regular number greater than or equal to target.
    Regular numbers are composites of the prime factors 2, 3, and 5.
    Also known as 5-smooth numbers or Hamming numbers, these are the optimal
    size for inputs to FFTPACK.
    Target must be a positive integer.
    """
    if target <= 6:
        return target

    # Quickly check if it's already a power of 2
    if not (target & (target-1)):
        return target

    match = float('inf')  # Anything found will be smaller
    p5 = 1
    while p5 < target:
        p35 = p5
        while p35 < target:
            # Ceiling integer division, avoiding conversion to float
            # (quotient = ceil(target / p35))
            quotient = -(-target // p35)

            # Quickly find next power of 2 >= quotient
            try:
                p2 = 2**((quotient - 1).bit_length())
            except AttributeError:
                # Fallback for Python <2.7
                p2 = 2**(len(bin(quotient - 1)) - 2)

            N = p2 * p35
            if N == target:
                return N
            elif N < match:
                match = N
            p35 *= 3
            if p35 == target:
                return p35
        if p35 < match:
            match = p35
        p5 *= 5
        if p5 == target:
            return p5
    if p5 < match:
        match = p5
    return match


def _fourier_transform_fft(t, ft, n='auto', axis=-1, conv=-1):
    r"""
    $f(\omega) = \int{dt} f(t) e^{-conv i\omega t}$
    times is assumed sorted and regularly spaced
    conv = 1 or -1 allows to change convention. Default is -1.
    """

    if (len(t) != ft.shape[axis]):
        raise ValueError, "coordinates should have the same length as values array on specified `axis`."
    if (conv not in [-1, +1]):
        raise ValueError, "`conv` must be -1 or +1."
    if n is None:
        n = len(t)
    elif n == 'auto':
        n = _next_regular(len(t))

    dt = t[1] - t[0]

    w = _fft.fftshift(_fft.fftfreq(n, dt))
    fw = _fft.fftshift(_fft.fft(ft, n=n, axis=axis), axes=axis)
    fw = _np.swapaxes(fw, -1, axis)

    w = 2 * _np.pi * conv * w[::conv]
    fw = fw[..., ::conv]
    fw[..., :] *= dt * _np.exp(-conv * 1j * w * t[0])

    return w, _np.swapaxes(fw, -1, axis)


def get_K_R_time(bin_times, kernel_binner):
    """
    Compute the retarded kernel K (a given series order) in time domain from a kernel binner.

    This does *not* multiply by i^n.
    This assumes the binner starts at t=0.
    """
    tmax = bin_times[0] + bin_times[-1]  # because binner starts at t=0
    times = bin_times - tmax
    bin_size = abs(times[1] - times[0])
    kernel = (kernel_binner[:, 0] + kernel_binner[:, 1]) / bin_size

    ### now we have the advanced kernel, transform it into the retarded one
    times = -times[::-1]
    kernel = _np.conj(kernel[::-1])

    return times, kernel


def get_G_R_omega(g0_R_omega, bin_times, kernel_binner, g0_is_vectorized=False):
    """
    Compute the retarded Green function (a given series order) in frequency domain from a kernel
    binner.

    FFT is used to move to the frequency domain.

    Arguments:
    - g0_R_omega: a function producing the unperturbed retarded Green function omega -> g_0^R(omega).
      This is e.g. G0ModelOmega(params).g0_dot_R.
    - bin_times: times of center of bins, obtained from the binner,
    - kernel_binner: the binner returned by keldy.
    - g0_is_vectorized (bool): tells if g0_R_omega is in vectorized form, e.g. using _np.vectorize.

    This does *not* multiply by i^n.
    This assumes the binner starts at t=0.

    Returns omegas, G_n^R(omega)
    """
    times, kernel = get_K_R_time(bin_times, kernel_binner)

    ### into frequency domain
    if not g0_is_vectorized:
        g0_R_omega = _np.vectorize(g0_R_omega)

    omegas, kernel = _fourier_transform_fft(times, kernel)
    gf = kernel * g0_R_omega(omegas)

    return omegas, gf


def get_G_R_omega_series(g0_R_omega, bin_times, kernel_binner_series, g0_is_vectorized=False):
    """
    Compute the retarded Green function series in frequency domain, from a series of kernel binners.

    This is a simple wrap around `get_G_R_omega` which is called for each order. In addition,
    order 0 is added at the beginning of the series.

    Returns omega and G^R_n(omega) as a (N,W)-array, N is number of orders, W number of frequencies.
    """
    gf_series = []

    if not g0_is_vectorized:
        g0_R_omega = _np.vectorize(g0_R_omega)

    for k, kernel_binner in enumerate(kernel_binner_series):
        omegas, gf = get_G_R_omega(g0_R_omega, bin_times, kernel_binner_series[k],
                                   g0_is_vectorized=True)
        gf_series.append(gf)

    gf_series.insert(0, g0_R_omega(omegas))
    return omegas, _np.array(gf_series, dtype=complex)


def get_G_R_omega_series_vect(g0_R_omega, bin_times, kernel_binner_series_list, g0_is_vectorized=False):
    """
    Compute a list of retarded Green function series in frequency domain, from a list of
    series of kernel binners.

    This is a simple wrap around `get_G_R_omega_series` which is called for each element of the list.

    Returns omega and G^R_n(omega) as a (M,N,W)-array, M is length of list, N number of orders,
    W number of frequencies.
    """

    gf_series_list = []

    if not g0_is_vectorized:
        g0_R_omega = _np.vectorize(g0_R_omega)

    for k, kernel_binner_series in enumerate(kernel_binner_series_list):
        omegas, gf_series = get_G_R_omega_series(g0_R_omega, bin_times, kernel_binner_series,
                                                 g0_is_vectorized=True)
        gf_series_list.append(gf_series)

    return omegas, _np.array(gf_series_list, dtype=complex)
