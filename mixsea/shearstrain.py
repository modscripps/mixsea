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


def strain_polynomial_fits(depth, t, SP, lon, lat, depth_bin, dz):
    """
    Calculate strain with a smooth :math:`N^2` profile from polynomial fits to windowed data.

    Parameters
    ----------
    depth : array-like
        CTD depth [m]
    t : array-like
        CTD in-situ temperature [ITS-90, degrees C]
    SP : array-like
        CTD practical salinity [psu]
    lon : array-like or float
        Longitude
    lat : array-like or float
        Latitude
    depth_bin : float
        Window centers
    dz : float
        Window size

    Returns
    -------
    strain : array-like
        Strain profile.
    depth_st : array-like
        Depth vector for `strain`.
    N2ref : array-like
        Smooth :math:`N^2` profile.
    """
    p = gsw.p_from_z(-depth, lat)
    SA = gsw.SA_from_SP(SP, p, lon, lat)
    CT = gsw.CT_from_t(SA, t, p)
    N2, Pbar = gsw.Nsquared(SA, CT, p, lat=lat)
    Zbar = -1 * gsw.z_from_p(Pbar, lat)
    isN = np.isfinite(N2)
    n2polyfit = np.zeros(N2.shape) * np.nan
    n2polyfit_mean = n2polyfit.copy()
    for iwin, zw in enumerate(depth_bin):
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
    depth_st = depth[:-1] + np.diff(depth[:2] / 2)
    # N2ref = n2polyfit_mean
    N2ref = interp1d(Zbar, n2polyfit_mean, bounds_error=False)(depth_st)
    strain = interp1d(Zbar, strain, bounds_error=False)(depth_st)
    return strain, depth_st, N2ref


def strain_adiabatic_leveling(depth, t, SP, lon, lat, bin_width):
    """
    Calculate strain with a smooth :math:`N^2` profile based on the adiabatic leveling method.

    Parameters
    ----------
    depth : array-like
        CTD depth [m]
    t : array-like
        CTD in-situ temperature [ITS-90, degrees C]
    SP : array-like
        CTD practical salinity [psu]
    lon : array-like or float
        Longitude
    lat : array-like or float
        Latitude

    Returns
    -------
    strain : array-like
        Strain profile.
    depth_st : array-like
        Depth vector for `strain`.
    N2ref : array-like
        Smooth :math:`N^2` profile.
    """
    p = gsw.p_from_z(-depth, lat)
    N2ref = nsq.adiabatic_leveling(
        p,
        SP,
        t,
        lon,
        lat,
        bin_width=bin_width,
        order=2,
        return_diagnostics=False,
        cap="both",
    )
    # N2ref = interp1d(z, N2ref)()
    SA = gsw.SA_from_SP(SP, p, lon, lat)
    CT = gsw.CT_from_t(SA, t, p)
    N2, Pbar = gsw.Nsquared(SA, CT, p, lat=lat)
    Zbar = -1 * gsw.z_from_p(Pbar, lat)
    N2ref = interp1d(p, N2ref)(Pbar)
    strain = (N2 - N2ref) / N2ref
    depth_st = depth[:-1] + np.diff(depth[:2] / 2)
    strain = interp1d(Zbar, strain, bounds_error=False)(depth_st)
    return strain, depth_st, N2ref


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


def nan_shearstrain(
    depth, t, SP, lon, lat, ladcp_u=None, ladcp_v=None, ladcp_depth=None, **kwargs
):
    """Compute krho and epsilon via shear/strain parameterization attempting to deal NaN values in the input data.

    See `shearstrain` for details.
    """
    depth = np.asarray(depth)
    t = np.asarray(t)
    SP = np.asarray(SP)

    # Look for NaNs in CTD.
    notnan = np.isfinite(SP) & np.isfinite(t) & np.isfinite(depth)
    isnan = ~notnan

    # Look for NaNs in shear/velocity. These do not determine the output shape,
    # therefore we can just get rid of the NaNs inside `kwargs`.
    if ladcp_u is not None:
        notnansh = (
            np.isfinite(ladcp_u) & np.isfinite(ladcp_v) & np.isfinite(ladcp_depth)
        )
        isnansh = ~notnansh
        if isnansh.sum() > 0:
            ladcp_u = ladcp_u[notnansh]
            ladcp_v = ladcp_v[notnansh]
            ladcp_depth = ladcp_depth[notnansh]

    if isnan.sum() == 0:  # If there are no NaNs in CTD data then return.
        return shearstrain(
            depth, t, SP, lon, lat, ladcp_u, ladcp_v, ladcp_depth, **kwargs
        )

    # Don't want to pass return_diagnostics twice.
    if "return_diagnostics" in kwargs:
        return_diagnostics = kwargs.pop("return_diagnostics")
    else:
        return_diagnostics = False

    eps_shst, krho_shst, diag = shearstrain(
        depth[notnan],
        t[notnan],
        SP[notnan],
        lon,
        lat,
        ladcp_u,
        ladcp_v,
        ladcp_depth,
        return_diagnostics=True,
        **kwargs,
    )

    if return_diagnostics:
        return eps_shst, krho_shst, diag
    else:
        return eps_shst, krho_shst


def shearstrain(
    depth,
    t,
    SP,
    lon,
    lat,
    ladcp_u=None,
    ladcp_v=None,
    ladcp_depth=None,
    m=None,
    depth_bin=None,
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
    depth : array-like
        CTD depth [m]
    t : array-like
        CTD in-situ temperature [ITS-90, degrees C]
    SP : array-like
        CTD practical salinity [psu]
    lon : array-like or float
        Longitude
    lat : array-like or float
        Latitude
    ladcp_u : array-like, optional
        LADCP velocity east-west component [m/s]. If not provided, only the
        strain solution with a fixed shear/strain ratio of 3 will be computed.
    ladcp_v : array-like, optional
        LADCP velocity north-south component [m/s]
    ladcp_depth : array-like, optional
        LADCP depth vector [m]
    m : array-like, optional
        Wavenumber vector to interpolate spectra onto
    depth_bin : array-like, optional
        Centers of windows over which spectra are computed. Defaults to
        np.arange(75, max(depth), 150). Note that windows are half-overlapping so
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
        ``"depth_bin"``
            Center points of depth windows (`array-like`).
        ``"Nmean"``
            Average N per depth window calculated as root-mean from N^2(`array-like`).

    Notes
    -----
    Adapted from Jen MacKinnon and Amy Waterhouse.
    """
    # Average lon, lat into one value if they are vectors.
    lon = np.nanmean(lon)
    lat = np.nanmean(lat)

    depth = np.asarray(depth)
    t = np.asarray(t)
    SP = np.asarray(SP)

    # Make sure there are no NaNs in the input data
    notnan = np.isfinite(SP) & np.isfinite(t) & np.isfinite(depth)
    isnan = ~notnan
    if not isnan.sum() == 0:
        raise ValueError(
            "No NaNs allowed in CTD data. Consider using `nan_shearstrain` instead."
        )

    # Can we work with velocity data?
    calcsh = False if ladcp_u is None else True

    # No NaNs in velocity/shear data
    if calcsh:
        ladcp_u = np.asarray(ladcp_u)
        ladcp_v = np.asarray(ladcp_v)
        ladcp_depth = np.asarray(ladcp_depth)
        notnansh = (
            np.isfinite(ladcp_u) & np.isfinite(ladcp_v) & np.isfinite(ladcp_depth)
        )
        isnansh = ~notnansh
        if not isnansh.sum() == 0:
            raise ValueError(
                "No NaNs allowed in LADCP data. Consider using `nan_shearstrain` instead."
            )

    # Coriolis parameter for this latitude
    f = np.absolute(gsw.f(lat))

    if calcsh:
        depth_sh = ladcp_depth

        # Calculate shear if necessary
        if ladcp_is_shear is False:
            uz = helpers.calc_shear(ladcp_u, depth_sh)
            vz = helpers.calc_shear(ladcp_v, depth_sh)
        else:
            uz = ladcp_u
            vz = ladcp_v

    # Create an evenly spaced wavenumber vector if none was provided
    if m is None:
        m = wavenumber_vector(w=300)

    # Generate depth bin vector if none provided
    if depth_bin is None:
        depth_bin = np.arange(75, np.max(depth), 150)
    else:
        # cut out any bins that won't hold any data
        depth_bin = depth_bin[depth_bin < np.max(depth)]
    nz = np.squeeze(depth_bin.shape)
    delz = np.mean(np.diff(depth_bin))

    # Calculate a smoothed N^2 profile and strain, either using 2nd order
    # polynomial fits to N^2 for each window (PF) or the adiabatic leveling
    # method (AL).
    if smooth == "PF":
        strain, depth_st, N2ref = strain_polynomial_fits(
            depth, t, SP, lon, lat, depth_bin, delz
        )
    elif smooth == "AL":
        strain, depth_st, N2ref = strain_adiabatic_leveling(
            depth, t, SP, lon, lat, bin_width=300
        )

    # Buoyancy-normalize shear
    if calcsh:
        N2ref_adcp = interp1d(depth_st, N2ref, bounds_error=False)(depth_sh)
        shear_un = uz / np.real(np.sqrt(N2ref_adcp))
        shear_vn = vz / np.real(np.sqrt(N2ref_adcp))
        shearn = shear_un + 1j * shear_vn

    # Remove nan in shear, strain
    if calcsh:
        iin = np.isfinite(shearn)
        shearn = shearn[iin]
        depth_sh = depth_sh[iin]
        N2ref_adcp = N2ref_adcp[iin]

    iin = np.isfinite(strain)
    strain = strain[iin]
    depth_st = depth_st[iin]

    N2 = N2ref[iin]

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
    Int_st = np.full(nz, np.nan)
    Int_sh = np.full(nz, np.nan)

    for iwin, zi in enumerate(depth_bin):
        zw = depth_bin[iwin]

        # Shear spectra
        if calcsh:
            iz = (depth_sh >= (zw - delz)) & (depth_sh <= (zw + delz))
            ig = ~np.isnan(shearn[iz])
            if ig.size > 10:
                dz = np.mean(np.diff(depth_sh))
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
        iz = (depth_st >= (zw - delz)) & (depth_st <= (zw + delz))
        ig = ~np.isnan(strain[iz])
        if ig.size > 10:
            dz = np.mean(np.diff(depth_st))
            _, _, Ptot, m0 = helpers.psd(strain[iz], dz, ffttype="t", detrend=True)
            # Compensation for first differencing
            H = np.sinc(m0 * dz / 2 / np.pi) ** 2
            Ptot = Ptot / H
            Ptot_st = interp1d(m0, Ptot, bounds_error=False)(m)
            P_strain[iwin, :] = Ptot_st
        else:
            P_strain[iwin, :] = np.nan
            Ptot_st = np.zeros_like(m) * np.nan

        # Mean stratification per segment
        Nm = np.sqrt(np.nanmean(N2[iz]))
        Nmean[iwin] = Nm

        # Shear cutoff wavenumber
        if calcsh:
            iimsh, Mmax_sh[iwin] = find_cutoff_wavenumber(
                Ptot_sh[m_include_sh], m[m_include_sh], sh_integration_limit
            )

        # Strain cutoff wavenumber
        iimst, Mmax_st[iwin] = find_cutoff_wavenumber(
            Ptot_st[m_include_st], m[m_include_st], st_integration_limit
        )

        if calcsh:
            # Integrate shear spectrum to obtain shear variance
            Ssh = np.trapz(Ptot_sh[iimsh], m[iimsh])
            Int_sh[iwin] = Ssh
            # GM shear variance
            Sshgm, Pshgm = gm_shear_variance(m, iimsh, Nm)

        # Integrate strain spectrum to obtain strain variance
        Sst = np.trapz(Ptot_st[iimst], m[iimst])
        Int_st[iwin] = Sst
        # GM strain variance
        Sstgm, Pstgm = gm_strain_variance(m, iimst, Nm)

        if calcsh:
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
            depth_bin=depth_bin,
            Nmean=Nmean,
            strain=strain,
            depth_st=depth_st,
            Int_sh=Int_sh,
            Int_st=Int_st,
        )
        return eps_shst, krho_shst, diag
    else:
        return eps_shst, krho_shst
