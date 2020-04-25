import gsw
import numpy as np


def eps_overturn(
    P, Z, T, S, lon, lat, dnoise=5e-4, alpha_sq=0.9, background_eps=np.nan
):
    """
    Calculate turbulent dissipation based on the Thorpe scale method.

    Parameters
    ----------
    P : array-like
        Pressure [dbar]
    Z : array-like
        Depth [m]
    T : array-like
        In-situ temperature [ITS90, °C]
    S : array-like
        Salinity [PSU]
    lon : float
        Longitude of observation
    lat : float
        Latitude of observation
    dnoise : float
        Noise level of density [kg/m^3] (default 5e-4)
    alpha_sq : float
        Square of proportionality constant between Thorpe and Ozmidov scale.
        (default 0.9)
    background_eps : float
        Background epsilon where no overturn detected. Defaults to np.nan.

    Returns
    -------
    Lt : array-like
      Thorpe length scale [m]
    eps : array-like
      Turbulent dissipation [W/kg]
    k : array-like
      Turbulent diffusivity [m^2/s]
    n2 : array-like
      Background stratification of each overturn detected [s^-2]
    dtdz : array-like
      Temperature gradient of each overturn [°/m]

    """
    # average lon, lat into one value if they are vectors
    lon = np.nanmean(lon)
    lat = np.nanmean(lat)

    # avoid error due to nan's in conditional statements
    # np.seterr(invalid="ignore")

    z0 = Z.copy()
    z0 = z0.astype("float")

    # populate output dict
    out = {}
    out["Lt"] = np.zeros_like(z0) * np.nan
    out["eps"] = np.zeros_like(z0) * np.nan
    out["k"] = np.zeros_like(z0) * np.nan
    out["n2"] = np.zeros_like(z0) * np.nan
    out["dtdz"] = np.zeros_like(z0) * np.nan

    # We need to calculate potential density and temperature at reference depths
    # If the profile is shallower than 1200 m, use only one depth range.
    # Otherwise, use several depth reference levels.
    if (np.nanmax(P) - np.nanmin(P)) > 1200:
        dref = 1000
        refd = np.arange(np.nanmin(P) + dref / 2, np.nanmax(P), dref)
    else:
        refd = (np.nanmin(P) + np.nanmax(P)) / 2
        dref = (np.nanmax(P) - np.nanmin(P)) + 1

    for refdi in refd:
        # Find non-NaNs
        x = np.squeeze(np.where(np.isfinite(T + S + P)))

        # Extract variables without the NaNs
        p = P[x].copy()
        z = Z[x].copy()
        z = z.astype("float")
        t = T[x].copy()
        s = S[x].copy()
        # cn2   = ctdn['n2'][x].copy()

        SA = gsw.SA_from_SP(s, t, lon, lat)
        CT = gsw.CT_from_t(SA, t, p)
        # PT = gsw.pt0_from_t(SA, t, p)

        # Calculate potential density
        sg = gsw.pot_rho_t_exact(SA, t, p, p_ref=refdi)

        # Create intermediate density profile
        D0 = sg[0]
        sgt = D0 - sg[0]
        n = sgt / dnoise
        n = np.fix(n)
        sgi = [D0 + n * dnoise]  # first element
        for i in np.arange(1, len(sg), 1):
            sgt = sg[i] - sgi[i - 1]
            n = sgt / dnoise
            n = np.fix(n)
            sgi.append(sgi[i - 1] + n * dnoise)
        sgi = np.array(sgi)

        # Sort (important to use mergesort here)
        # Ds = np.sort(sgi, kind="mergesort")
        Is = np.argsort(sgi, kind="mergesort")

        # Sort temperature profile as well for calculation of dT/dz
        # Thetas = np.sort(PT, kind='mergesort')
        # ThIs = np.argsort(PT, kind='mergesort')

        # Calculate Thorpe length scale
        TH = z[Is] - z
        cumTH = np.cumsum(TH)

        # make sure there are any overturns
        if np.sum(cumTH) > 2:

            aa = np.where(cumTH > 2)[0]
            blocks = _consec_blocks(aa, combine_gap=1)

            # Sort temperature and salinity based on the density sorting index
            # for calculating the buoyancy frequency
            # PTs = PT[Is]
            SAs = SA[Is]
            CTs = CT[Is]

            # Loop over detected overturns and calculate Thorpe Scales, N2
            # and dT/dz over the overturn region
            THsc = np.zeros_like(z) * np.nan
            N2 = np.zeros_like(z) * np.nan
            # CN2  = np.ones_like(z)*np.nan
            DTDZ = np.zeros_like(z) * np.nan

            for iostart, ioend in zip(blocks[:, 0], blocks[:, 1]):
                idx = np.arange(iostart, ioend + 1, 1)
                sc = np.sqrt(np.mean(np.square(TH[idx])))
                # ctdn2 = np.nanmean(cn2[idx])
                # Buoyancy frequency calculated over the overturn from sorted
                # profiles. Go beyond overturn.
                n2, Np = gsw.Nsquared(
                    SAs[[iostart - 1, ioend + 1]],
                    CTs[[iostart - 1, ioend + 1]],
                    p[[iostart - 1, ioend + 1]],
                    lat,
                )
                # Fill depth range of the overturn with the Thorpe scale
                THsc[idx] = sc
                # Fill depth range of the overturn with N^2
                N2[idx] = n2

                # Fill depth range of the overturn with local temperature gradient
                if iostart > 0:
                    PTov = CTs[iostart - 1 : ioend + 1]
                    zov = z[iostart - 1 : ioend + 1]
                else:
                    PTov = CTs[iostart : ioend + 1]
                    zov = z[iostart : ioend + 1]
                local_dtdz = (np.min(PTov) - np.max(PTov)) / (np.max(zov) - np.min(zov))
                DTDZ[idx] = local_dtdz

            # Find data for this reference depth range
            iz = np.squeeze(
                np.where(((p > refdi - dref / 2) & (p <= refdi + dref / 2)))
            )
            # There shouldn't be any negative N^2 the way we calculate it -
            # but only over the current depth range!
            # Also, exclude any nan's to avoid warnings.
            ni = np.where(~np.isnan(N2[iz]))[0]
            assert np.all(N2[iz[ni]] > 0)
            # Calculate epsilon
            THepsilon = np.zeros_like(N2) * np.nan
            THepsilon[iz[ni]] = alpha_sq * THsc[iz[ni]] ** 2 * np.sqrt(N2[iz[ni]]) ** 3
            # THepsilon[N2 <= 0] = np.nan
            THk = np.zeros_like(N2) * np.nan
            THk[iz[ni]] = 0.2 * THepsilon[iz[ni]] / N2[iz[ni]]

            # Pick only data for this reference depth range
            out["eps"][x[iz]] = THepsilon[iz]
            out["k"][x[iz]] = THk[iz]
            out["n2"][x[iz]] = N2[iz]
            out["Lt"][x[iz]] = THsc[iz]
            out["dtdz"][x[iz]] = DTDZ[iz]

        # Fill with background epsilon
        ni = np.where(np.isnan(out["eps"][x]))
        out["eps"][x[ni]] = background_eps

    return out["Lt"], out["eps"], out["k"], out["n2"], out["dtdz"]


def _consec_blocks(idx=None, combine_gap=0, combine_run=0):
    """
    block_idx = consec_blocks(idx,combine_gap=0, combine_run=0)

    Routine that returns the start and end indexes of the consecutive blocks
    of the index array (idx). The second argument combines consecutive blocks
    together that are separated by <= combine. This is useful when you want
    to perform some action on the n number of data points either side of a
    gap, say, and don't want that action to be effected by a neighbouring
    gap.

    From Glenn Carter, University of Hawaii
    """
    if idx.size == 0:
        return np.array([])

    # Sort the index data and remove any identical points
    idx = np.unique(idx)

    # Find the block boundaries
    didx = np.diff(idx)
    ii = np.concatenate(((didx > 1).nonzero()[0], np.atleast_1d(idx.shape[0] - 1)))

    # Create the block_idx array
    block_idx = np.zeros((ii.shape[0], 2), dtype=int)
    block_idx[0, :] = [idx[0], idx[ii[0]]]
    for c in range(1, ii.shape[0]):
        block_idx[c, 0] = idx[ii[c - 1] + 1]
        block_idx[c, 1] = idx[ii[c]]

    # Find the gap between and combine blocks that are closer together than
    # the combine_gap threshold
    gap = (block_idx[1:, 0] - block_idx[0:-1, 1]) - 1
    if np.any(gap <= combine_gap):
        count = 0
        new_block = np.zeros(block_idx.shape, dtype=int)
        new_block[0, 0] = block_idx[0, 0]
        for ido in range(block_idx.shape[0] - 1):
            if gap[ido] > combine_gap:
                new_block[count, 1] = block_idx[ido, 1]
                count += 1
                new_block[count, 0] = block_idx[ido + 1, 0]
        new_block[count, 1] = block_idx[-1, 1]
        block_idx = new_block[: count + 1, :]

    # Combine any runs that are shorter than the combine_run threshold
    runlength = block_idx[:, 1] - block_idx[:, 0]
    if np.any(runlength <= combine_run):
        count = 0
        new_block = np.zeros(block_idx.shape, dtype=int)
        for ido in range(block_idx.shape[0]):
            if runlength[ido] > combine_run:
                new_block[count, :] = block_idx[ido, :]
                count += 1
        block_idx = new_block[:count, :]

    return np.atleast_2d(block_idx)
