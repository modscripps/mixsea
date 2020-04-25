import numpy as np
import pytest
from scipy import signal

from mixsea import helpers


def test_psd_output_dims():
    """Make sure output dimensions of psd work as expected."""

    # generate input matrix
    m = np.arange(0, 30)
    M = m.size
    n = np.linspace(0, 12 * np.pi, 1000)
    N = n.size
    g = np.tile(np.sin(n), (M, 1))
    mu, sigma = 0, 0.1
    noise = np.random.normal(mu, sigma, (M, N))
    g = g + noise * 10 + n / 10
    g = signal.detrend(g)

    # also a row input vector
    g1 = g[0, :]

    for ffttype in ["p", "s", "c", "t"]:
        for gin in [g1, g1[1:], g1.transpose(), g1[1:].transpose()]:
            Pcw, Pccw, Ptot, omega = helpers.psd(gin, dx=1, axis=-1, ffttype=ffttype)
            assert omega.shape == Ptot.shape
        for gin in [g, g[:, 1:], g.transpose(), g[:, 1:].transpose()]:
            Pcw, Pccw, Ptot, omega = helpers.psd(gin, dx=1, axis=-1, ffttype=ffttype)
            assert omega.shape[0] == Ptot.shape[1]
            Pcw, Pccw, Ptot, omega = helpers.psd(gin, dx=1, axis=0, ffttype=ffttype)
            assert omega.shape[0] == Ptot.shape[1]

    # test something different than the standard window
    Pcw, Pccw, Ptot, omega = helpers.psd(g, dx=1, axis=-1, ffttype="t", window="hann")
    assert omega.shape[0] == Ptot.shape[1]

    # test window size
    Pcw, Pccw, Ptot, omega = helpers.psd(g, dx=1, axis=-1, ffttype="t", tser_window=100)
    assert omega.shape[0] == Ptot.shape[1]

    # overlap larger than window should raise
    with pytest.raises(ValueError):
        Pcw, Pccw, Ptot, omega = helpers.psd(
            g, dx=1, axis=-1, ffttype="t", tser_window=100, tser_overlap=200
        )
