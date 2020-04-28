import gsw
import numpy as np
from scipy.interpolate import interp1d

from . import helpers


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
    m_include_sh=None,
    m_include_st=None,
    ladcp_is_shear=False,
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
        Wavenumber integration range for shear spectra. Array must consist of indices
        or boolans to index m. Defaults to first 4 wavenumbers.
    m_include_st : array-like, optional
        Wavenumber integration range for strain spectra. Array must consist of indices
        or boolans to index m. Defaults to first 4 wavenumbers.
    ladcp_is_shear : bool, optional
        Indicate whether LADCP data is velocity or shear.
        Defaults to False (velocity).

    Returns
    -------
    P_shear : array-like
        Matrix of shear spectra for each depth window
    P_strain : array-like
        Matrix of shear spectra for each depth window
    Mmax_sh : array-like
        Cutoff wavenumber (kc)
    Mmax_st : array-like
        Cutoff wavenubmer used for strain only calculation
    Rwtot : array-like
        Shear/strain ratio used, computed from spectra unless specificed in input
    krho_shst : array-like
        krho calculated from both shear and strain spectra
    krho_st : array-like
        krho calculated from strain only
    eps_shst : array-like
        Epsilon calculated from both shear and strain spectra
    eps_st : array-like
        Epsilon calculated from strain only
    m : array-like
        Wavenumber vector
    z_bin : array-like
        Center points of depth windows
    Nmean : array-like
        Average N per depth window calculated as root-mean from n2 above
    """
    # average lon, lat into one value if they are vectors
    lon = np.nanmean(lon)
    lat = np.nanmean(lat)

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

    if z_bin is None:
        z_bin = np.arange(75, np.max(z), 150)
    else:
        # cut out any bins that won't hold any data
        z_bin = z_bin[z_bin < np.max(z)]
    nz = np.squeeze(z_bin.shape)

    # Calculate buoyancy frequency by fitting a 2nd order polynomial
    SA = gsw.SA_from_SP(sa, pr, lon, lat)
    CT = gsw.CT_from_t(SA, te, pr)
    N2, Pbar = gsw.Nsquared(SA, CT, pr, lat=lat)
    # Convert output pressure to depth
    Zbar = -1 * gsw.z_from_p(Pbar, lat)
    n2polyfit = np.zeros(N2.shape) * np.nan
    n2polyfit_mean = n2polyfit.copy()
    for jj in np.array(range(nz - 1)):
        if jj < nz - 2:
            ij = np.squeeze(np.where(((Zbar > z_bin[jj]) & (Zbar <= z_bin[jj + 2]))))
        elif jj == nz - 1:
            ij = np.squeeze(np.where(((Zbar > z_bin[jj]) & (Zbar <= z_bin[jj + 1]))))
        p = np.polyfit(Zbar[ij], N2[ij], deg=2)
        pv = np.polyval(p, Zbar[ij])
        n2polyfit_mean[ij] = np.tile(np.mean(pv), ij.shape)
        n2polyfit[ij] = pv

    # Extraplolate over NaN's at the bottom
    n2polyfit_mean = helpers.extrapolate_data(n2polyfit_mean)
    n2polyfit = helpers.extrapolate_data(n2polyfit)

    # Interpolate to LADCP depths
    n2polyfit_mean_adcp = interp1d(Zbar, n2polyfit_mean, bounds_error=False)(dd)
    # Interpolate to CTD depths
    n2polyfit_mean_ctd = interp1d(Zbar, n2polyfit_mean, bounds_error=False)(z)

    # Buoyancy-normalize shear
    shear_un = uz / np.real(np.sqrt(n2polyfit_mean_adcp))
    shear_vn = vz / np.real(np.sqrt(n2polyfit_mean_adcp))
    shearn = shear_un + 1j * shear_vn
    z_sh = dd.copy()

    # Calculate strain
    strain = (N2 - n2polyfit) / n2polyfit_mean
    z_st = z[:-1] + np.diff(z[:2] / 2)

    # Create an evenly spaced wavenumber vector
    if m is None:
        m = np.arange(2 * np.pi / 300, 2 * np.pi / 10, 2 * np.pi / 600)

    # Remove nan in shear, strain
    iin = np.isfinite(shearn)
    shearn = shearn[iin]
    z_sh = z_sh[iin]
    n2polyfit_mean_adcp = n2polyfit_mean_adcp[iin]

    iin = np.isfinite(strain)
    strain = strain[iin]
    z_st = z_st[iin]

    # Run shear/strain calculation
    (
        P_shear,
        P_strain,
        Mmax_sh,
        Mmax_st,
        Rwtot,
        krho_shst,
        krho_st,
        Nmean,
    ) = compute_shearstrain_krho(
        shearn,
        z_sh,
        strain,
        z_st,
        n2polyfit_mean_adcp,
        lat,
        m,
        z_bin,
        m_include_sh=m_include_sh,
        m_include_st=m_include_st,
    )

    # Calculate turbulent dissipation from krho
    n2 = interp1d(pr, n2polyfit_mean_ctd)(z_bin)
    eps_shst = 1 / 0.2 * krho_shst * n2
    eps_st = 1 / 0.2 * krho_st * n2

    return (
        P_shear,
        P_strain,
        Mmax_sh,
        Mmax_st,
        Rwtot,
        krho_shst,
        krho_st,
        eps_shst,
        eps_st,
        m,
        z_bin,
        Nmean,
    )


def compute_shearstrain_krho(
    shearn,
    z_sh,
    strain,
    z_st,
    n2,
    lat,
    m,
    z_bin,
    Rwavg=None,
    m_include_sh=None,
    m_include_st=None,
):
    """
    Compute vertical diffusivity based on the shear/strain parameterization.

    Adapted from Jen MacKinnon.

    Parameters
    ----------
    shearn : array-like
        Vertical shear profile (du/dz+i*dv/dz) normalized by a SMOOTHED version of N
    z_sh : array-like
        Depth corresponding to shearn
    strain : array-like
        Strain profile
    z_st : array-like
        Depth corresponding to strain
    n2 : array-like
        SMOOTHED version of N^2 interpolated onto SHEAR depth vector
    lat : array-like or float
        Profile latitude
    m : array-like
        Wavenumber vector to interpolate spectra onto
    z_bin : array-like
        Centers of windows over which spectra are computed, typically some-
        thing like 50:150:max(z), where 50 is the base of the mixed layer.
        Note that windows are half-overlapping so the 150 spacing above means
        each window is 300 m tall.
    Rwavg : array-like, optional
        Shear/strain ratio of z_bin.shape if you'd like to specify it (based
        for example on a station average), otherwise computed using shear and
        strain spectra for each spectral window.
    m_include_sh : array-like, optional
        Wavenumber integration range for shear spectra. Array must consist of indices
        or boolans to index m. Defaults to first 4 wavenumbers.
    m_include_st : array-like, optional
        Wavenumber integration range for strain spectra. Array must consist of indices
        or boolans to index m. Defaults to first 4 wavenumbers.

    Returns
    -------
    P_shear : array-like
        Matrix of shear spectra for each depth window
    P_strain : array-like
        Matrix of shear spectra for each depth window
    Mmax_sh : array-like
        Cutoff wavenumber (kc)
    Mmax_st : array-like
        Cutoff wavenubmer used for strain only calculation
    Rwtot : array-like
        Shear/strain ratio used, computed from spectra unless specificed in input
    krho_shst : array-like
        krho calculated from both shear and strain spectra
    krho_st : array-like
        krho calculated from strain only
    Nmean : array-like
        Average N per depth window calculated as root-mean from n2 above
    """
    # Convert wavenumber includes in case they are given as range().
    if m_include_sh is not None:
        m_include_sh = np.asarray(m_include_sh)
    if m_include_st is not None:
        m_include_st = np.asarray(m_include_st)

    delz = np.mean(np.diff(z_bin))
    nz = z_bin.size

    K0 = 0.05 * 1e-4
    N0 = 5.24e-3

    P_shear = np.zeros((nz, m.size)) * np.nan
    P_strain = P_shear.copy()
    Mmax_sh = np.zeros((nz)) * np.nan
    Rwtot = Mmax_sh.copy()
    krho_shst = Mmax_sh.copy()
    krho_st = krho_shst.copy()
    Mmax_st = Mmax_sh.copy()
    Nmean = Mmax_sh.copy()

    f = np.absolute(gsw.f(lat))
    f30 = gsw.f(30)

    for iwin, zi in enumerate(z_bin):
        zw = z_bin[iwin]
        iz = (z_sh >= (zw - delz)) & (z_sh <= (zw + delz))
        nn = np.sqrt(np.nanmean(n2[iz]))
        Nmean[iwin] = nn
        Jf = f * np.arccosh(nn / f) / f30 / np.arccosh(N0 / f30)

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

        # Find cutoff wavenumber based on SHEAR spectra.
        # See methods in Gregg et al 2003 for factor 0.66
        if m_include_sh is not None:
            iim = m_include_sh
        else:
            # print('selecting first four wavenumbers for strain integration')
            iim = np.array(range(4))
        specsum = np.cumsum(Ptot_sh * np.mean(np.diff(m)))
        iim2 = np.where(np.less(specsum[iim], 0.66, where=np.isfinite(specsum[iim])))[0]
        if iim2.size > 0:
            iim = iim[iim2]
        if iim.size > 1:
            Mmax_sh[iwin] = np.max(m[iim])
        else:
            Mmax_sh[iwin] = np.nan

        # Shear/strain ratio
        if Rwavg:
            Rw = Rwavg[iwin]
        else:
            if np.any(np.isfinite(Ptot_sh[iim])) & np.any(np.isfinite(Ptot_st[iim])):
                Rw = np.nanmean(Ptot_sh[iim]) / np.nanmean(Ptot_st[iim])
                Rw = 1.01 if Rw < 1.01 else Rw
            else:
                Rw = np.nan
        Rwtot[iwin] = Rw

        hRw = 3 * (Rw + 1) / (2 * np.sqrt(2) * Rw * np.sqrt(Rw - 1))

        # Gm shear spectra
        Pgm = (
            (3 * np.pi * 6.3e-5 * 1300 * 3 / 2)
            * m ** 2
            / (m + 3 * np.pi / 1300 * nn / N0) ** 2
        )

        krho_shst[iwin] = (
            K0 * np.sum(Ptot_sh[iim]) ** 2 / np.sum(Pgm[iim]) ** 2 * hRw * Jf
        )

        # find cutoff wavenumber based on STRAIN spectra
        if m_include_st is not None:
            iim = m_include_st
        else:
            # print('selecting first four wavenumbers for strain integration')
            iim = np.array(range(4))
        specsum = np.cumsum(Ptot_st[iim] * np.mean(np.diff(m[iim])))
        iim2 = np.where(np.less(specsum, 0.22, where=np.isfinite(specsum)))[0]
        iim = iim[iim2]
        if iim.size > 1:
            # Do not go to wavenumbers corresponding to lambda_z < 5m
            if np.max(m[iim] / 2 / np.pi) > 0.2:
                iim = np.where(m / 2 / np.pi < 0.2)[0]
            Mmax_st[iwin] = np.max(m[iim])
            # Gm shear spectra - note shear/strain ratio 3 in here when comparing to shear
            Pgm_st = (
                (np.pi * 6.3e-5 * 1300 * 3 / 2)
                * m ** 2
                / (m + 3 * np.pi / 1300 * nn / N0) ** 2
            )
            # Use assumed shear/strain ratio of 3
            Rw = 3
            # for strain, different for shear:
            h2Rw = 1 / 6 / np.sqrt(2) * Rw * (Rw + 1) / np.sqrt(Rw - 1)
            krho_st[iwin] = (
                K0 * np.sum(Ptot_st[iim]) ** 2 / np.sum(Pgm_st[iim]) ** 2 * h2Rw * Jf
            )
        else:
            Mmax_st[iwin] = np.nan
            krho_st[iwin] = np.nan

    return P_shear, P_strain, Mmax_sh, Mmax_st, Rwtot, krho_shst, krho_st, Nmean
