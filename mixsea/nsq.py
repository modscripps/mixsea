import gsw
import numpy as np
from scipy import interpolate, signal


def nsqfcn(s, t, p, p0, dp, lon, lat):
    """Calculate square of buoyancy frequency [rad/s]^2 for profile of
    temperature, salinity and pressure.

    The Function: (1) low-pass filters t,s,p over dp,
                  (2) linearly interpolates  t and s onto pressures, p0,
                      p0+dp, p0+2dp, ....,
                  (3) computes upper and lower potential densities at
                      p0+dp/2, p0+3dp/2,...,
                  (4) converts differences in potential density into nsq
                  (5) returns NaNs if the filtered pressure is not
                      monotonic.

    If you want to have the buoyancy frequency in [cyc/s] then
    calculate sqrt(n2)./(2*pi). For the period in [s] do sqrt(n2).*2.*pi

    Adapted from Gregg and Alford.

    Gunnar Voet
    gvoet@ucsd.edu

    Parameters
    ----------
    s : float
        Salinity
    t : float
        In-situ temperature
    p : float
        Pressures
    p0 : float
        Lower bound (start) pressure for output values (not important...)
    dp : float
        Pressure interval of output data
    lon : float
        Longitude of observation
    lat : float
        Latitude of observation

    Returns
    -------
    n2 : Buoyancy frequency squared in (rad/s)^2
    pout : Pressure vector for n2

    """
    G = 9.80655

    # Make sure data has dtype np.ndarray
    if type(s) is not np.ndarray:
        s = np.array(s)
    if type(p) is not np.ndarray:
        p = np.array(p)
    if type(t) is not np.ndarray:
        t = np.array(t)

    # Delete negative pressures
    xi = np.where(p >= 0)
    p = p[xi]
    s = s[xi]
    t = t[xi]

    # Exclude nan in t and s
    xi = np.where((~np.isnan(s)) & (~np.isnan(t)))
    p = p[xi]
    s = s[xi]
    t = t[xi]

    # Put out all nan if no good data left
    if ~p.any():
        n2 = np.nan
        pout = np.nan

    # Reverse order of upward profiles
    if p[-1] < p[0]:
        p = p[::-1]
        t = t[::-1]
        s = s[::-1]

    # Low pass filter temp and salinity to match specified dp
    dp_data = np.diff(p)
    dp_med = np.median(dp_data)
    # [b,a]=butter(4,2*dp_med/dp); %causing problems...
    a = 1
    b = np.hanning(2 * np.floor(dp / dp_med))
    b = b / np.sum(b)

    tlp = signal.filtfilt(b, a, t)
    slp = signal.filtfilt(b, a, s)
    plp = signal.filtfilt(b, a, p)

    # Check that p is monotonic
    if np.all(np.diff(plp) >= 0):
        pmin = plp[0]
        pmax = plp[-1]

        # # Sort density if opted for
        #   if sort_dens
        #     rho = sw_pden(slp,tlp,plp,plp);
        #     [rhos, si] = sort(rho,'ascend');
        #     tlp = tlp(si);
        #     slp = slp(si);
        #   end

        while p0 <= pmin:
            p0 = p0 + dp

        # End points of nsq window
        pwin = np.arange(p0, pmax, dp)
        ft = interpolate.interp1d(plp, tlp)
        t_ep = ft(pwin)
        fs = interpolate.interp1d(plp, slp)
        s_ep = fs(pwin)
        # Determine the number of output points
        (npts,) = t_ep.shape

        # Compute pressures at center points
        pout = np.arange(p0 + dp / 2, np.max(pwin), dp)

        # Compute potential density of upper window pts at output pressures
        sa_u = gsw.SA_from_SP(s_ep[0:-1], t_ep[0:-1], lon, lat)
        pd_u = gsw.pot_rho_t_exact(sa_u, t_ep[0:-1], pwin[0:-1], pout)

        # Compute potential density of lower window pts at output pressures
        sa_l = gsw.SA_from_SP(s_ep[1:], t_ep[1:], lon, lat)
        pd_l = gsw.pot_rho_t_exact(sa_l, t_ep[1:], pwin[1:], pout)

        # Compute buoyancy frequency squared
        n2 = G * (pd_l - pd_u) / (dp * pd_u)

    else:
        print("  filtered pressure not monotonic")
        n2 = np.nan
        pout = np.nan

    return n2, pout


def adiabatic_leveling(
    P,
    S,
    T,
    lon,
    lat,
    bin_width=100.0,
    order=1,
    return_diagnostics=False,
    cap=None,
):
    """Generate smooth buoyancy frequency profile by adiabatic leveling.

    Parameters
    ----------
    P : 1-D ndarray
        Pressure [dbar]
    S : 1-D ndarray
        Practical salinity [-]
    T : 1-D ndarray
        Temperature [degrees C]
    lon : float
        Longitude [-180...+360]
    lat : float
        Latitude [-90...+90]
    bin_width : float, optional
        Pressure bin width [dbar]
    order : int, optional
        Degree of polynomial fit.
    return_diagnostics : bool, optional
        Flag to return additional argument pcoefs. False by default.
    cap : {None, 'left', 'right', 'both'}, optional
        Flag to change proceedure at ends of array where bins may be partially
        filled. None by default, meaning they are included. Can also specify
        'left', 'right' or 'both' to cap method before partial bins.

    Returns
    -------
    N2_ref : 1-D ndarray
        Reference buoyancy frequency [s-2]
    pcoefs : 2-D ndarray
        Fitting coefficients, returned only when the flag return_diagnostics is set
        True.

    Notes
    -----
    Based on method of Bray and Fofonoff (1981) :cite:`Bray1981`.

    """
    valid = np.isfinite(P) & np.isfinite(S) & np.isfinite(T)
    valid = np.squeeze(np.argwhere(valid))
    P_, S_, T_ = P[valid], S[valid], T[valid]

    flip = False
    if (np.diff(P_) < 0).all():
        flip = True
        P_ = np.flipud(P_)
        S_ = np.flipud(S_)
        T_ = np.flipud(T_)
    elif (np.diff(P_) < 0).any():
        raise ValueError("P must be monotonically increasing/decreasing.")

    i1 = np.searchsorted(P_, P_ - bin_width / 2.0)
    i2 = np.searchsorted(P_, P_ + bin_width / 2.0)

    if cap is None:
        Nd = P_.size
    elif cap == "both":
        icapl = i2[0]
        icapr = i1[-1]
    elif cap == "left":
        icapl = i2[0]
        icapr = i2[-1]
    elif cap == "right":
        icapl = i1[0]
        icapr = i1[-1]
    else:
        raise ValueError(
            "The argument cap must be either None, 'both', 'left'" " or 'right'"
        )

    if cap is not None:
        i1 = i1[icapl:icapr]
        i2 = i2[icapl:icapr]
        valid = valid[icapl:icapr]
        Nd = icapr - icapl

    dimax = np.max(i2 - i1)

    Pb = np.full((dimax, Nd), np.nan)
    Sb = np.full((dimax, Nd), np.nan)
    Tb = np.full((dimax, Nd), np.nan)

    for i in range(Nd):
        imax = i2[i] - i1[i]
        Pb[:imax, i] = P_[i1[i] : i2[i]]
        Sb[:imax, i] = S_[i1[i] : i2[i]]
        Tb[:imax, i] = T_[i1[i] : i2[i]]

    Pbar = np.nanmean(Pb, axis=0)

    SAb = gsw.SA_from_SP(Sb, Pb, lon, lat)

    rho = gsw.pot_rho_t_exact(SAb, Tb, Pb, Pbar)
    sv = 1.0 / rho

    rhobar = np.nanmean(rho, axis=0)

    p = np.full((order + 1, Nd), np.nan)

    for i in range(Nd):
        imax = i2[i] - i1[i]
        p[:, i] = np.polyfit(Pb[:imax, i], sv[:imax, i], order)

    g = gsw.grav(lat, Pbar)
    # The factor 1e-4 is needed for conversion from dbar to Pa.
    if order == 1:
        N2 = -1e-4 * rhobar ** 2 * g ** 2 * p[0, :]
    elif order == 2:
        N2 = -1e-4 * rhobar ** 2 * g ** 2 * (p[1, :] + 2 * Pbar * p[0, :])
    elif order == 3:
        N2 = (
            -1e-4
            * rhobar ** 2
            * g ** 2
            * (p[2, :] + 2 * Pbar * p[1, :] + 3 * Pbar ** 2 * p[0, :])
        )
    else:
        raise ValueError("Fits are only included up to 3rd order.")

    N2_ref = np.full_like(P, np.nan)
    pcoef = np.full((order + 1, P.size), np.nan)
    if flip:
        N2_ref[valid] = np.flipud(N2)
        pcoef[:, valid] = np.fliplr(p)
    else:
        N2_ref[valid] = N2
        pcoef[:, valid] = p

    if return_diagnostics:
        return N2_ref, pcoef
    else:
        return N2_ref
