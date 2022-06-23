import pathlib

import numpy as np
from scipy import fftpack, signal
from scipy.interpolate import interp1d


def denan(data):
    ii = np.squeeze(np.argwhere(~np.isnan(data)))
    return data[ii], ii


def calc_shear(u, depth):
    """Calculate vertical shear from velocity data.

    Parameters
    ----------
    u : array-like
        Velocity data.
    depth : array-like
        Depth vector.

    Returns
    -------
    array-like
        Vertical shear at depths of `depth`.
    """
    f = interp1d(
        depth[0:-1] + np.diff(depth[:2]) / 2,
        np.diff(u) / np.diff(depth),
        bounds_error=False,
    )
    return f(depth)


def extrapolate_data(data):
    """Extrapolate over NaNs near the bottom."""
    out = data.copy()
    tmp = np.squeeze(np.where(~np.isnan(data)))
    max_tmp = np.max(tmp)
    if np.isfinite(max_tmp) & np.isnan(data[-1]):
        maxind = max_tmp
        out[maxind + 1 :] = data[maxind]
    return out


def psd(
    g,
    dx,
    axis=1,
    ffttype="p",
    detrend=True,
    window="hamming",
    tser_window=None,
    tser_overlap=None,
):
    """
    Compute power spectral density.

    Adapted from Jen MacKinnon.

    Parameters
    ----------
    g : array-like
        Real or complex input data of size [M * 1], [1 * M] or [M * N]
    dx : float
        Distance or time between entries in g
    axis : int, optional
        Axis along which to fft: 0 = columns, 1 = rows (default)
    ffttype : str, optional
        Flag for fft type

        - 'p' for periodic boundary conditions = pure fft (default)

        - 'c' or 's' if the input is a sum of cosines or sines only in which
          case the input is extended evenly or oddly before ffting.

        - 't' for time series, or any other non-exactly periodic series in
          which the data should be windowed and filtered before computing the
          periodogram. In this case you may also specify:

          * tser_window: an integer that gives the length of each window
            the series should be broken up into for filtering and
            computing the periodogram.  Default is length(g).

          * tser_overlap: an integer equal to the lengh of points
            overlap between sucessive windows. Default = tser_window/2,
            which means 50% overlap.
    detrend : bool, optional
        Detrend along dim by removing a linear fit using scipy.detrend().
        Defaults to True.
    window : str, optional
        Window type. Default 'hamming'. See scipy.signal.get_window() for
        window options.

    Returns
    -------
    Pcw : array-like
        Clockwise power spectral density
    Pccw : array-like
        Counter-clockwise power spectral density
    Ptot : array-like
        Total = ccw + cw psd one-sided spectra
    omega: array-like
        Vector of wavenumbers. For plotting use frequency f = omega/(2*pi).
    """
    M0 = g.shape
    g = np.atleast_2d(g)
    if detrend:
        g = signal.detrend(g)

    # FFT and welch act on rows by default. If we want to calculate column-wise,
    # simply transpose the input matrix.
    if axis == 0 & len(M0) > 1:
        g = np.transpose(g)
    M = g.shape

    # If sin or cos transform needed, appropriately extend time series
    if ffttype == "s":
        g_ext = np.concatenate(
            (g, -1 * np.flip(g[:, 1 : int(M[1]) - 1], axis=1)), axis=1
        )
        g = g_ext
        M = g.shape
    if ffttype == "c":
        g_ext = np.concatenate((g, np.flip(g[:, 1 : int(M[1]) - 1], axis=1)), axis=1)
        g = g_ext
        M = g.shape

    # Setup frequency vectors
    df = 1 / (M[1] * dx)
    domega = 2 * np.pi * df
    # full frequency output length M
    f_full = np.linspace(0, (M[1] - 1) * df, num=M[1], endpoint=True)
    if np.remainder(M[1], 2) == 0:  # even -> length M/2+1
        f, step = np.linspace(
            0, M[1] * df / 2, num=int(M[1] / 2 + 1), endpoint=True, retstep=True
        )
    else:  # odd -> length (M+1)/2
        f, step = np.linspace(
            0, (M[1] - 1) * df / 2, num=int((M[1] + 1) / 2), endpoint=True, retstep=True
        )
    assert step - df < df / 1000
    omega = 2 * np.pi * f

    # Compute power spectra using fft
    if ffttype in ["p", "c", "s"]:
        P0 = fftpack.fft(g, axis=-1)
        # Normalize by wavenumber in RADIANS
        Pxx = P0 * np.conj(P0) / M[1] / (M[1] * domega)

    if ffttype == "t":
        if tser_window is None:
            # default one segment
            tser_window = np.array(np.fix(M[1]))
        if tser_overlap is None:
            # default 50% overlap
            tser_overlap = np.fix(tser_window / 2)

        f0, Pxx0 = signal.welch(
            g,
            fs=1 / dx,
            axis=-1,
            nperseg=tser_window,
            noverlap=tser_overlap,
            window=window,
            return_onesided=False,
        )

        # One side of f0 is shifted to negative values - shift everything to
        # positive.
        f0 = np.fft.fftshift(f0) + np.absolute(f0).max()

        # Interpolate to f_full frequency vector
        Pxx0 = interp1d(f0, Pxx0, bounds_error=False, axis=-1)(f_full)
        Pxx = Pxx0.copy()
        Pxx = Pxx / 2 / np.pi  # normalize by radial wavenumber/frequency

    # Separate into cw and ccw spectra
    if np.remainder(M[1], 2) == 0:
        # even lengths, divide 0 and nyquist freq between ccw and cw
        Pccw = np.concatenate(
            (
                np.reshape(0.5 * Pxx[:, 0], (-1, 1)),
                Pxx[:, 1 : int((M[1]) / 2)],
                np.reshape(0.5 * Pxx[:, int(M[1] / 2) + 1], (-1, 1)),
            ),
            axis=1,
        )
        Pcw = np.concatenate(
            (
                np.reshape(0.5 * Pxx[:, 0], (-1, 1)),
                np.flip(Pxx[:, int(M[1] / 2 + 1) : M[1]], axis=1),
                np.reshape(0.5 * Pxx[:, int(M[1] / 2) + 2], (-1, 1)),
            ),
            axis=1,
        )
        Ptot = Pccw + Pcw
    else:
        # odd lengths, divide 0 freq between ccw and cw
        Pccw = np.concatenate(
            (np.reshape(0.5 * Pxx[:, 0], (-1, 1)), Pxx[:, 1 : int((M[1] + 1) / 2)]),
            axis=1,
        )
        Pcw = np.concatenate(
            (
                np.reshape(0.5 * Pxx[:, 0], (-1, 1)),
                np.flip(Pxx[:, int(((M[1] + 3) / 2)) - 1 : M[1]], axis=1),
            ),
            axis=1,
        )
        Ptot = Pccw + Pcw

    Ptot = np.squeeze(Ptot)
    Pcw = np.squeeze(Pcw)
    Pccw = np.squeeze(Pccw)

    # transpose back if axis=0
    if axis == 0 & len(M0) > 1:
        Ptot = Ptot.transpose()
        Pcw = Pcw.transpose()
        Pccw = Pccw.transpose()

    return Pcw, Pccw, Ptot, omega


def read_testfile(filename):
    """
    Read ascii file with test data

    The ascii file is expected to have one header line with variable names,
    preceded by a #. The variable names are used as keys in the output dict.

    Parameters
    ----------
    filename : str or pathlib.Path Path to testfile

    Returns
    -------
    dict Dictionary with test data
    """
    data = np.loadtxt(filename, comments="#", delimiter=",")
    # read variable names from comment line and divide data into a dict
    with open(filename) as f:
        vars = f.readline().strip()
        v = [x.strip() for x in vars[2:].split(",")]
    return {vi: data[:, i] for i, vi in enumerate(v)}


def read_ctd_testfile():
    p = pathlib.Path(__file__).resolve()
    ctdfile = p.parent / "tests" / "data" / "ctd_test_data.csv"
    ctd = read_testfile(ctdfile)
    return ctd


def read_ladcp_testfile():
    p = pathlib.Path(__file__).resolve()
    ladcpfile = p.parent / "tests" / "data" / "ladcp_test_data.csv"
    ladcp = read_testfile(ladcpfile)
    return ladcp
