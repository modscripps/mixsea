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
