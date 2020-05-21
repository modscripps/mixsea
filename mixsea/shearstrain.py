import gsw
import numpy as np
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d

from . import helpers, nsq


def wavenumber_vector(w):
    r"""
    Generate wavenumber vector

    Wavenumber vector runs from
    :math:`\frac{2 \pi}{\textrm{w}}` to
    :math:`\frac{2 \pi}{10}` in increments of
    :math:`\frac{2 \pi}{\textrm{w}}`.

    Parameters
    ----------
    w : float
        window size [m]

    Returns
    -------
    m : array
        Wavenumber vector
    """
    return np.arange(2 * np.pi / w, 2 * np.pi / 10, 2 * np.pi / w)


def strain_polynomial_fits(s, t, p, z, lon, lat, zbin, dz):
    """
    Calculate strain with a smooth N^2 profile from polynomial fits to windowed data.

    Parameters
    ----------
    s : array-like
        CTD salinity [psu]
    t : array-like
        CTD in-situ temperature [ITS-90, degrees C]
    p : array-like
        CTD pressure [dbar]
    z : array-like
        CTD depth [m]
    lat : array-like or float
        Latitude
    lon : array-like or float
        Longitude
    zbin : float
        Window centers
    dz : float
        Window size

    Returns
    -------
    strain : array-like
        Strain profile.
    z_st : array-like
        Depth vector for `strain`.
    N2ref : array-like
        Smooth N^2 profile.
    """
    SA = gsw.SA_from_SP(s, p, lon, lat)
    CT = gsw.CT_from_t(SA, t, p)
    N2, Pbar = gsw.Nsquared(SA, CT, p, lat=lat)
    Zbar = -1 * gsw.z_from_p(Pbar, lat)
    isN = np.isfinite(N2)
    n2polyfit = np.zeros(N2.shape) * np.nan
    n2polyfit_mean = n2polyfit.copy()
    for iwin, zw in enumerate(zbin):
        ij = (Zbar >= (zw - dz)) & (Zbar < (zw + dz))
        ij2 = (Zbar >= (zw - dz / 2)) & (Zbar < (zw + dz / 2))
        ij = ij * isN
        p = np.polyfit(Zbar[ij], N2[ij], deg=2)
        pv = np.polyval(p, Zbar[ij2])
        n2polyfit_mean[ij2] = np.mean(pv)
        n2polyfit[ij2] = pv
    # Extraplolate over NaN's at the bottom
    n2polyfit_mean = helpers.extrapolate_data(n2polyfit_mean)
    n2polyfit = helpers.extrapolate_data(n2polyfit)
    # Calculate strain
    strain = (N2 - n2polyfit) / n2polyfit_mean
    z_st = z[:-1] + np.diff(z[:2] / 2)
    # N2ref = n2polyfit_mean
    N2ref = interp1d(Zbar, n2polyfit_mean, bounds_error=False)(z_st)
    strain = interp1d(Zbar, strain, bounds_error=False)(z_st)
    return strain, z_st, N2ref


def strain_adiabatic_leveling(s, t, p, z, lon, lat, bin_width):
    """
    Calculate strain with a smooth N^2 profile based on the adiabatic leveling method.

    Parameters
    ----------
    s : array-like
        CTD salinity [psu]
    t : array-like
        CTD in-situ temperature [ITS-90, degrees C]
    p : array-like
        CTD pressure [dbar]
    z : array-like
        CTD depth [m]
    lat : array-like or float
        Latitude
    lon : array-like or float
        Longitude

    Returns
    -------
    strain : array-like
        Strain profile.
    z_st : array-like
        Depth vector for `strain`.
    N2ref : array-like
        Smooth N^2 profile.
    """
    N2ref = nsq.adiabatic_leveling(
        p,
        s,
        t,
        lon,
        lat,
        bin_width=bin_width,
        order=2,
        return_diagnostics=False,
        cap="both",
    )
    # N2ref = interp1d(z, N2ref)()
    SA = gsw.SA_from_SP(s, p, lon, lat)
    CT = gsw.CT_from_t(SA, t, p)
    N2, Pbar = gsw.Nsquared(SA, CT, p, lat=lat)
    Zbar = -1 * gsw.z_from_p(Pbar, lat)
    N2ref = interp1d(p, N2ref)(Pbar)
    strain = (N2 - N2ref) / N2ref
    z_st = z[:-1] + np.diff(z[:2] / 2)
    strain = interp1d(Zbar, strain, bounds_error=False)(z_st)
    return strain, z_st, N2ref


def find_cutoff_wavenumber(P, m, integration_limit, lambda_min=5):
    """
    Find cutoff wavenumber for spectrum integration.

    Parameters
    ----------
    P : array-like
        Spectrum
    m : array-like
        Vertical wavenumber vector
    integration_limit : float
        Variance limit for integration
    lambda_min : float, optional
        Minimum vertical wavelength lambda_z=2pi/m, by default 5

    Returns
    -------
    iim : array-like
        Integration range as indexer to `P`.
    """
    specsum = cumtrapz(P, m)
    specsum = np.insert(specsum, 0, 0)
    iim = np.flatnonzero(
        np.less(specsum, integration_limit, where=np.isfinite(specsum))
    )
    if iim.size == 1:  # need at least two data points for integration
        iim = np.append(iim, 1)
    if iim.size > 1:
        # Do not go to wavenumbers corresponding to lambda_z < lambda_min meter
        if np.max(m[iim] / 2 / np.pi) > 0.2:
            iim = np.where(m / 2 / np.pi < 0.2)[0]
        return iim, np.max(m[iim])
    else:
        return iim, np.nan


def latitude_correction(f, N):
    r"""
    Latitudinal correction term

    Parameters
    ----------
    f : float
        Coriolis parameter
    N : float
        Buoyancy frequency

    Returns
    -------
    L : float
        Latitudinal correction

    Notes
    -----
    Calculates the latitudinal dependence term as described in Gregg et al.
    (2003) :cite:`Gregg2003`:

    .. math::

        L(\theta, N) = \frac{f \cosh^{-1}(N/f)}{f_{30^{\circ}} \cosh^{-1}(N_0/f_{30^\circ})}

    with Coriolis parameter at 30° latitude :math:`f_{30^\circ}` and reference
    GM buoyancy frequency :math:`N_0=5.24\times10^{-3}\,\mathrm{s}^{-1}`.
    """
    # Coriolis parameter at 30 degrees latitude
    f30 = gsw.f(30)  # rad s-1
    # GM model reference stratification:
    N0 = 5.24e-3  # rad s-1
    f = np.abs(f)
    return f * np.arccosh(N / f) / (f30 * np.arccosh(N0 / f30))


def gm_shear_variance(m, iim, N):
    r"""
    GM model shear variance

    Parameters
    ----------
    m : array-like
        Vertical wavenumber vector [rad/m]
    iim : array-like
        Wavenumber integration range, indexer to m
    N : float
        Local buoyancy frequency [s^-1]

    Returns
    -------
    Sgm : float
        GM shear variance normalized by N^2 [1/m^2].
    Pgm : array-like
        GM shear spectrum for wavenumber range `m`.

    Notes
    -----
    Returns GM shear variance normalized by buoyancy frequency by integrating
    `Pgm` over a wavenumber range of the GM shear spectrum  as presented in
    Kunze et al. (2006) :cite:`Kunze2006` eq. 6:

    .. math::

        \frac{_{GM}\left< V_z^2\right>}{\overline{N}^2} = \frac{3 \pi E_0 b j_\ast}{2} \int_{m_\mathrm{min}}^{m_\mathrm{max}} \frac{m^2 dm}{(m + m_\ast)^2}

    with :math:`j_\ast=3`, :math:`E_0=6.3\times10^{-5}`, :math:`N_0=5.2\times
    10^{-3}` rad/s, :math:`b=1300` m, and

    .. math::

        m_\ast = \frac{\overline{N}}{N_0}\frac{\pi j_\ast}{b}

    See also
    --------
    gm_strain_variance : GM strain variance
    """
    N0 = 5.24e-3  # reference buoyancy frequency = 3 cph
    b = 1300  # thermocline scale depth
    jstar = 3
    E0 = 6.3e-5  # GM energy level
    Pgm = (
        (3 * np.pi * E0 * b * jstar / 2)
        * m ** 2
        / (m + jstar * np.pi / b * N / N0) ** 2
    )
    # integrate
    Sgm = np.trapz(y=Pgm[iim], x=m[iim])
    return Sgm, Pgm


def gm_strain_variance(m, iim, N):
    r"""
    GM model strain variance.

    Parameters
    ----------
    m : array-like
        Vertical wavenumber vector [rad/m].
    iim : array-like
        Wavenumber integration range, indexer to `m`.
    N : float
        Local buoyancy frequency [s^-1].

    Returns
    -------
    Sgm : float
        GM strain variance normalized by N^2 [1/m^2].
    Pgm : array-like
        GM strain spectrum for wavenumber range `m`.

    Notes
    -----
    Returns GM strain variance by integrating `Pgm` over a wavenumber range of
    the GM shear spectrum  as presented in Kunze et al. (2006)
    :cite:`Kunze2006` eq. 10:

    .. math::

        \frac{_{GM}\left< V_z^2\right>}{\overline{N}^2} = \frac{\pi E_0 b j_\ast}{2} \int_{m_\mathrm{min}}^{m_\mathrm{max}} \frac{m^2 dm}{(m + m_\ast)^2}

    with :math:`j_\ast=3`, :math:`E_0=6.3\times10^{-5}`, :math:`N_0=5.2\times
    10^{-3}` rad/s, :math:`b=1300` m, and

    .. math::

        m_\ast = \frac{\overline{N}}{N_0}\frac{\pi j_\ast}{b}

    Note that this corresponds to the buoyancy-normalized GM shear variance
    divided by 3.

    See also
    --------
    gm_shear_variance : GM shear variance

    """
    N0 = 5.24e-3  # reference buoyancy frequency = 3 cph
    b = 1300  # thermocline scale depth
    jstar = 3
    E0 = 6.3e-5  # GM energy level
    Pgm = (np.pi * E0 * b * jstar / 2) * m ** 2 / (m + jstar * np.pi / b * N / N0) ** 2
    # integrate
    Sgm = np.trapz(y=Pgm[iim], x=m[iim])
    return Sgm, Pgm


def shearstrain(
    s,
    t,
    p,
    z,
    lat,
    lon,
    ladcp_u,
    ladcp_v,
    ladcp_z,
    m=None,
    z_bin=None,
    m_include_sh=np.arange(4),
    m_include_st=np.arange(4),
    ladcp_is_shear=False,
    smooth="AL",
    sh_integration_limit=0.66,
    st_integration_limit=0.22,
    return_diagnostics=False,
):
    """
    Compute krho and epsilon from CTD/LADCP data via the shear/strain parameterization.

    Parameters
    ----------
    s : array-like
        CTD salinity [psu]
    t : array-like
        CTD in-situ temperature [ITS-90, degrees C]
    p : array-like
        CTD pressure [dbar]
    z : array-like
        CTD depth [m]
    lat : array-like or float
        Latitude
    lon : array-like or float
        Longitude
    ladcp_u : array-like
        LADCP velocity east-west component [m/s]
    ladcp_v : array-like
        LADCP velocity north-south component [m/s]
    ladcp_z : array-like
        LADCP depth vector [m]
    m : array-like, optional
        Wavenumber vector to interpolate spectra onto
    z_bin : array-like, optional
        Centers of windows over which spectra are computed. Defaults to
        np.arange(75, max(z), 150). Note that windows are half-overlapping so
        the 150 spacing above means each window is 300 m tall.
    m_include_sh : array-like, optional
        Wavenumber integration range for shear spectra. Array must consist of
        indices or boolans to index m. Defaults to first 4 wavenumbers.
    m_include_st : array-like, optional
        Wavenumber integration range for strain spectra. Array must consist of
        indices or boolans to index m. Defaults to first 4 wavenumbers.
    ladcp_is_shear : bool, optional
        Indicate whether LADCP data is velocity or shear.
        Defaults to False (velocity).
    smooth : {'AL', 'PF'}, optional
        Select type of N^2 smoothing and subsequent strain calculation.
        Defaults to adiabatic leveling.
    sh_integration_limit : float
        Shear variance level for determining integration cutoff wavenumber.
        Defaults to 0.66, compare Gregg et al. (2003).
    st_integration_limit : float
        Strain variance level for determining integration cutoff wavenumber.
        Defaults to 0.22, compare Gregg et al. (2003).
    return_diagnostics : bool, optional
        Default is False. If True, this function will return a dictionary
        containing variables such as shear spectra, shear/strain ratios,

    Returns
    -------
    eps_shst : array-like
        Epsilon calculated from both shear and strain spectra
    krho_shst : array-like
        krho calculated from both shear and strain spectra
    diag : dict, optional
        Dictionary of diagnostic variables, set return with the
        `return_diagnostics' argument. `diag` holds the following variables:

        ``"P_shear"``
            Matrix of shear spectra for each depth window (`array-like`).
        ``"P_strain"``
            Matrix of strain spectra for each depth window
        ``"Mmax_sh"``
            Cutoff wavenumber kc (`array-like`).
        ``"Mmax_st"``
            Cutoff wavenubmer used for strain only calculation (`array-like`).
        ``"Rwtot"``
            Shear/strain ratio used, computed from spectra unless specificed in
            input (`array-like`).
        ``"krho_st"``
            krho calculated from strain only (`array-like`).
        ``"eps_st"``
            Epsilon calculated from strain only (`array-like`).
        ``"m"``
            Wavenumber vector (`array-like`).
        ``"z_bin"``
            Center points of depth windows (`array-like`).
        ``"Nmean"``
            Average N per depth window calculated as root-mean from N^2(`array-like`).

    Notes
    -----
    Adapted from Jen MacKinnon and Amy Waterhouse.
    """
    # average lon, lat into one value if they are vectors
    lon = np.nanmean(lon)
    lat = np.nanmean(lat)

    # Coriolis parameter for this latitude
    f = np.absolute(gsw.f(lat))

    sa, ii = helpers.denan(s)
    te = t[ii]
    pr = p[ii]
    z = z[ii]

    u, ii = helpers.denan(ladcp_u)
    v = ladcp_v[ii]
    dd = ladcp_z[ii]

    # Calculate shear
    if ladcp_is_shear is False:
        uz = helpers.calc_shear(u, dd)
        vz = helpers.calc_shear(v, dd)
    else:
        uz = u
        vz = v

    # Create an evenly spaced wavenumber vector if none was provided
    if m is None:
        m = wavenumber_vector(w=300)

    # Generate depth bin vector if none provided
    if z_bin is None:
        z_bin = np.arange(75, np.max(z), 150)
    else:
        # cut out any bins that won't hold any data
        z_bin = z_bin[z_bin < np.max(z)]
    nz = np.squeeze(z_bin.shape)
    delz = np.mean(np.diff(z_bin))

    # Calculate a smoothed N^2 profile and strain, either using 2nd order
    # polynomial fits to N^2 for each window (PF) or the adiabatic leveling
    # method (AL).
    if smooth == "PF":
        strain, z_st, N2ref = strain_polynomial_fits(
            sa, te, pr, z, lon, lat, z_bin, delz
        )
    elif smooth == "AL":
        strain, z_st, N2ref = strain_adiabatic_leveling(
            sa, te, pr, z, lon, lat, bin_width=300
        )

    # Interpolate N2ref to LADCP depths
    N2ref_adcp = interp1d(z_st, N2ref, bounds_error=False)(dd)

    # Buoyancy-normalize shear
    shear_un = uz / np.real(np.sqrt(N2ref_adcp))
    shear_vn = vz / np.real(np.sqrt(N2ref_adcp))
    shearn = shear_un + 1j * shear_vn
    z_sh = dd.copy()

    # Remove nan in shear, strain
    iin = np.isfinite(shearn)
    shearn = shearn[iin]
    z_sh = z_sh[iin]
    N2ref_adcp = N2ref_adcp[iin]

    iin = np.isfinite(strain)
    strain = strain[iin]
    z_st = z_st[iin]

    N2 = N2ref_adcp

    # Convert wavenumber integration range to np.array in case it is something
    # else:
    m_include_sh = np.asarray(m_include_sh)
    m_include_st = np.asarray(m_include_st)

    N0 = 5.24e-3  # (3 cph)
    K0 = 0.05 * 1e-4
    eps0 = 7.8e-10  # Waterman et al. 2014
    # eps0 = 7.9e-10 # Polzin et al. 1995
    # eps0 = 6.73e-10 # Gregg et al. 2003

    P_shear = np.full((nz, m.size), np.nan)
    P_strain = P_shear.copy()
    Mmax_sh = np.full(nz, np.nan)
    Mmax_st = Mmax_sh.copy()
    Rwtot = np.full(nz, np.nan)
    krho_shst = np.full(nz, np.nan)
    krho_st = np.full(nz, np.nan)
    eps_shst = np.full(nz, np.nan)
    eps_st = np.full(nz, np.nan)
    Nmean = np.full(nz, np.nan)

    for iwin, zi in enumerate(z_bin):
        zw = z_bin[iwin]
        iz = (z_sh >= (zw - delz)) & (z_sh <= (zw + delz))
        Nm = np.sqrt(np.nanmean(N2[iz]))
        Nmean[iwin] = Nm

        # Shear spectra
        ig = ~np.isnan(shearn[iz])
        if ig.size > 10:
            dz = np.mean(np.diff(z_sh))
            _, _, Ptot, m0 = helpers.psd(shearn[iz], dz, ffttype="t", detrend=True)
            # Compensation for first differencing
            H = np.sinc(m0 * dz / 2 / np.pi) ** 2
            Ptot = Ptot / H
            Ptot_sh = interp1d(m0, Ptot, bounds_error=False)(m)
            P_shear[iwin, :] = Ptot_sh
        else:
            P_shear[iwin, :] = np.nan
            Ptot_sh = np.zeros_like(m) * np.nan

        # Strain spectra
        iz = (z_st >= (zw - delz)) & (z_st <= (zw + delz))
        ig = ~np.isnan(strain[iz])
        if ig.size > 10:
            dz = np.mean(np.diff(z_st))
            _, _, Ptot, m0 = helpers.psd(strain[iz], dz, ffttype="t", detrend=True)
            # Compensation for first differencing
            H = np.sinc(m0 * dz / 2 / np.pi) ** 2
            Ptot = Ptot / H
            Ptot_st = interp1d(m0, Ptot, bounds_error=False)(m)
            P_strain[iwin, :] = Ptot_st
        else:
            P_strain[iwin, :] = np.nan
            Ptot_st = np.zeros_like(m) * np.nan

        # Shear cutoff wavenumber
        iimsh, Mmax_sh[iwin] = find_cutoff_wavenumber(
            Ptot_sh[m_include_sh], m[m_include_sh], sh_integration_limit
        )

        # Strain cutoff wavenumber
        iimst, Mmax_st[iwin] = find_cutoff_wavenumber(
            Ptot_st[m_include_st], m[m_include_st], st_integration_limit
        )

        # Integrate shear spectrum to obtain shear variance
        Ssh = np.trapz(Ptot_sh[iimsh], m[iimsh])
        # GM shear variance
        Sshgm, Pshgm = gm_shear_variance(m, iimsh, Nm)

        # Integrate strain spectrum to obtain strain variance
        Sst = np.trapz(Ptot_st[iimst], m[iimst])
        # GM strain variance
        Sstgm, Pstgm = gm_strain_variance(m, iimst, Nm)

        # Shear/strain ratio normalized by GM. Factor 3 corrects for the ratio
        # of GM shear to strain = 3 N^2.
        Rw = 3 * (Ssh / Sshgm) / (Sst / Sstgm)
        Rwtot[iwin] = Rw
        # Avoid negative square roots in hRw below
        Rw = 1.01 if Rw < 1.01 else Rw

        # Shear/strain parameterization
        hRw = 3 * (Rw + 1) / (2 * np.sqrt(2) * Rw * np.sqrt(Rw - 1))
        krho_shst[iwin] = (
            K0 * (Ssh ** 2 / Sshgm ** 2) * (hRw * latitude_correction(f, Nm))
        )
        eps_shst[iwin] = (
            eps0
            * (Nm ** 2 / N0 ** 2)
            * (Ssh ** 2 / Sshgm ** 2)
            * (hRw * latitude_correction(f, Nm))
        )

        # Strain only parameterization
        # Use assumed shear/strain ratio of 3
        Rw = 3
        h2Rw = 1 / 6 / np.sqrt(2) * Rw * (Rw + 1) / np.sqrt(Rw - 1)
        krho_st[iwin] = (
            K0 * (Sst ** 2 / Sstgm ** 2) * (h2Rw * latitude_correction(f, Nm))
        )
        eps_st[iwin] = (
            eps0
            * (Nm ** 2 / N0 ** 2)
            * (Sst ** 2 / Sstgm ** 2)
            * (h2Rw * latitude_correction(f, Nm))
        )

    if return_diagnostics:
        diag = dict(
            eps_st=eps_st,
            krho_st=krho_st,
            P_shear=P_shear,
            P_strain=P_strain,
            Mmax_sh=Mmax_sh,
            Mmax_st=Mmax_st,
            Rwtot=Rwtot,
            m=m,
            z_bin=z_bin,
            Nmean=Nmean,
        )
        return eps_shst, krho_shst, diag
    else:
        return eps_shst, krho_shst
