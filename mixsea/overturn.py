import gsw
import numpy as np


def nan_eps_overturn(
    P,
    T,
    S,
    lon,
    lat,
    dnoise=5e-4,
    alpha=0.95,
    Roc=0.2,
    background_eps=np.nan,
    use_ip=True,
    Nsq_method="teosp1",
    overturns_from_CT=False,
    return_diagnostics=False,
):
    """
    Calculate turbulent dissipation based on the Thorpe scale method attempting to deal NaN values in the input data.

    Parameters
    ----------
    P : array-like
        Pressure [dbar]
    T : array-like
        In-situ temperature [ITS90, °C]
    S : array-like
        Salinity [PSU]. Can be a single constant value. This may be convenient if only temperature data
        are available.
    lon : float
        Longitude of observation
    lat : float
        Latitude of observation
    dnoise : float, optional
        Noise level of density [kg/m^3]. Default is 5e-4.
    alpha : float, optional
        Constant of proportionality between Thorpe and Ozmidov scale. Default is 0.95.
    Roc : float, optional
        Critical value of the overturn ratio Ro. An overturn will be considered
        noise if Ro < Roc.
    background_eps : float, optional
        Background epsilon where no overturn detected. Defaults to numpy.nan.
    use_ip : bool, optional
        Sets whether to use the intermediate profile method. Default is True. If True,
        the dnoise parameter is passed as the `accuracy' argument of the intermediate
        profile method.
    Nsq_method : string, optional
        Method for calculation of buoyancy frequency. Default is 'teosp1'. Options are 'bulk',
        'endpt', 'teos' and 'teosp1'.
    overturns_from_CT : bool, optional
        If true, overturning patches will be diagnosed from the conservative temperature CT, 
        instead of potential density. Default is False.
    return_diagnostics : bool, optional
        Default is False. If True, this function will return a dictionary containing
        variables such as the Thorpe scale Lt, etc.

    Returns
    -------
    eps : array-like
        Turbulent dissipation [W/kg]
    n2 : array-like
        Background stratification of each overturn detected [s^-2]
    diag : dict, optional
        Dictionary of diagnositc variables, set return with the `return_diagnostics' argument.


    """
    P = np.asarray(P)
    T = np.asarray(T)
    S = np.asarray(S)

    if S.size == 1:
        S = np.full_like(P, S)

    # Find non-NaNs
    notnan = np.isfinite(P) & np.isfinite(T) & np.isfinite(S)
    isnan = ~notnan

    if isnan.sum() == 0:  # If there are no NaNs then return.
        return eps_overturn(
            P,
            T,
            S,
            lon,
            lat,
            dnoise=dnoise,
            alpha=alpha,
            Roc=Roc,
            background_eps=background_eps,
            use_ip=use_ip,
            Nsq_method=Nsq_method,
            overturns_from_CT=overturns_from_CT,
            return_diagnostics=return_diagnostics,
        )

    eps = np.full_like(P, np.nan)
    n2 = np.full_like(P, np.nan)

    eps[notnan], n2[notnan], diag = eps_overturn(
        P[notnan],
        T[notnan],
        S[notnan],
        lon,
        lat,
        dnoise=dnoise,
        alpha=alpha,
        Roc=Roc,
        background_eps=background_eps,
        use_ip=use_ip,
        Nsq_method=Nsq_method,
        overturns_from_CT=overturns_from_CT,
        return_diagnostics=True,
    )

    if return_diagnostics:
        # Replace nans in diagnostics if the size and shape seems right:
        Nnotnans = notnan.sum()
        for key in diag:
            if (np.size(diag[key]) == Nnotnans) & (np.ndim(diag[key]) == 1):
                ar = np.full_like(P, np.nan)
                ar[notnan] = diag[key]
                diag[key] = ar  # This will wipe out the old item.

        return eps, n2, diag

    else:
        return eps, n2


def eps_overturn(
    P,
    T,
    S,
    lon,
    lat,
    dnoise=5e-4,
    alpha=0.95,
    Roc=0.2,
    background_eps=np.nan,
    use_ip=True,
    Nsq_method="teosp1",
    overturns_from_CT=False,
    return_diagnostics=False,
):
    """
    Calculate turbulent dissipation based on the Thorpe scale method. This function cannot handle
    NaNs in the input data, but there is another called `nan_eps_overturn' that attempts to.

    Parameters
    ----------
    P : array-like
        Pressure [dbar]
    T : array-like
        In-situ temperature [ITS90, °C]
    S : float or array-like
        Salinity [PSU]. Can be a single constant value. This may be convenient if only temperature data
        are available.
    lon : float
        Longitude of observation
    lat : float
        Latitude of observation
    dnoise : float, optional
        Noise level of density [kg/m^3]. Default is 5e-4.
    alpha : float, optional
        Constant of proportionality between Thorpe and Ozmidov scale. Default is 0.95.
    Roc : float, optional
        Critical value of the overturn ratio Ro. An overturn will be considered
        noise if Ro < Roc.
    background_eps : float, optional
        Background epsilon where no overturn detected. Defaults to numpy.nan.
    use_ip : bool, optional
        Sets whether to use the intermediate profile method. Default is True. If True,
        the dnoise parameter is passed as the `accuracy' argument of the intermediate
        profile method.
    Nsq_method : string, optional
        Method for calculation of buoyancy frequency. Default is 'teosp1'. Options are 'bulk',
        'endpt', 'teos' and 'teosp1'.
    overturns_from_CT : bool, optional
        If true, overturning patches will be diagnosed from the conservative temperature CT, 
        instead of potential density. Default is False.
    return_diagnostics : bool, optional
        Default is False. If True, this function will return a dictionary containing
        variables such as the Thorpe scale Lt, etc.

    Returns
    -------
    eps : array-like
        Turbulent dissipation [W/kg]
    n2 : array-like
        Background stratification of each overturn detected [s^-2]
    diag : dict, optional
        Dictionary of diagnositc variables, set return with the `return_diagnostics' argument.


    """
    P = np.asarray(P)
    T = np.asarray(T)
    S = np.asarray(S)

    if not np.all(np.isclose(np.maximum.accumulate(P), P)):
        raise ValueError(
            "It appears that depth is not monotonically increasing, please fix."
        )

    if S.size == 1:
        S = np.full_like(P, S)

    if not (P.size == T.size == S.size):
        raise ValueError(
            "Input array sizes do not match. P.size = {}, T.size = {}, S.size = {}".format(
                P.size, T.size, S.size
            )
        )

    if not any(s == Nsq_method for s in ["teosp1", "teos", "endpt", "bulk"]):
        raise ValueError(
            "The 'Nsq_method' argument must be 'teosp1', 'teos', 'endpt' or 'bulk'."
        )

    ndata = P.size

    # Estimate depth from pressure.
    depth = -gsw.z_from_p(P, lat)
    dz = 0.5 * (depth[2:] - depth[:-2])  # 'width' of each data point
    dz = np.hstack((dz[0], dz, dz[-1]))  # assume width of first and last data point

    SA = gsw.SA_from_SP(S, T, lon, lat)
    CT = gsw.CT_from_t(SA, T, P)

    # Initialise arrays for diagnostic variables and flags.
    diag = {}
    diagvar = [
        "eps",
        "n2",
        "Lt",
        "thorpe_disp",
        "dens",
        "dens_sorted",
        "Ro",
        "SA_sorted",
        "CT_sorted",
    ]
    for var in diagvar:
        diag[var] = np.full_like(P, np.nan)

    if use_ip and not overturns_from_CT:
        diag["dens_ip"] = np.full_like(P, np.nan)

    if use_ip and overturns_from_CT:
        diag["CT_ip"] = np.full_like(P, np.nan)

    flagvar = ["noise_flag", "n2_flag", "ends_flag", "Ro_flag"]
    for var in flagvar:
        diag[var] = np.full_like(P, False, dtype=bool)

    # Potential density is only meaningful near the reference pressure. For a deep profile
    # we may need to select several reference pressures. To do so, we find the pressure
    # bins that best contain the data
    pbinwidth = 1000.0  # In future we could have this as an argument.
    pbinmin = np.floor(P.min() / pbinwidth) * pbinwidth
    pbinmax = np.ceil(P.max() / pbinwidth) * pbinwidth
    pbins = np.arange(pbinmin, pbinmax + pbinwidth, pbinwidth)
    p_refs = 0.5 * (
        pbins[1:] + pbins[:-1]
    )  # Use mid point pressure as reference pressure.
    nbins = p_refs.size

    # Loop over pressure bins.
    for idx_bin in range(nbins):

        dens = gsw.pot_rho_t_exact(SA, T, P, p_ref=p_refs[idx_bin])

        if overturns_from_CT:
            # Temperature normally decreases towards the bottom which would mean the
            # find_overturns algorithm thinks the whole water column is unstable! Minus fixes that.
            q = -CT
        else:
            q = dens

        if use_ip:  # Create intermediate density profile
            q_ip = intermediate_profile(
                q, acc=dnoise, hinge=1000, kind="down"
            )  # TODO: make hinge optional
            sidx, patches = find_overturns(
                q_ip, combine_gap=0
            )  # Also, combine gap should be a parameter...
        else:
            sidx, patches = find_overturns(q, combine_gap=0)

        # If there are no overturns, move on to the next pressure bin.
        if not np.any(patches):
            continue

        # Thorpe displacements
        thorpe_disp = depth[sidx] - depth

        q_sorted = q[sidx]

        # Sort other quantities based on the sorting indices.
        dens_sorted = dens[sidx]
        SA_sorted = SA[sidx]
        CT_sorted = CT[sidx]

        # Temporary arrays.
        Lt = np.full_like(P, np.nan)
        n2 = np.full_like(P, np.nan)
        Ro = np.full_like(P, np.nan)
        noise_flag = np.full_like(P, False, dtype=bool)
        n2_flag = np.full_like(P, False, dtype=bool)
        ends_flag = np.full_like(P, False, dtype=bool)
        Ro_flag = np.full_like(P, False, dtype=bool)

        for patch in patches:
            # Get patch indices.
            i0 = patch[0]
            i1 = patch[1] + 1  # Need +1 for last point in overturn.
            pidx = np.arange(i0, i1 + 1, 1)

            # Thorpe scale is the root mean square thorpe displacement.
            Lto = np.sqrt(np.mean(np.square(thorpe_disp[pidx])))
            Lt[pidx] = Lto

            # Estimate the buoyancy frequency.
            if Nsq_method == "teos":
                n2o, _ = gsw.Nsquared(
                    SA_sorted[[i0, i1]], CT_sorted[[i0, i1]], P[[i0, i1]], lat,
                )
            elif Nsq_method == "teosp1":
                # Go beyond overturn. Need to add 1 for this, unless end or beginning.
                addi = 0 if i1 == ndata - 1 else 1
                subi = 0 if i0 == 0 else 1

                n2o, _ = gsw.Nsquared(
                    SA_sorted[[i0 - subi, i1 + addi]],
                    CT_sorted[[i0 - subi, i1 + addi]],
                    P[[i0 - subi, i1 + addi]],
                    lat,
                )
            elif Nsq_method == "bulk":
                g = gsw.grav(lat, P[pidx].mean())
                densanom = dens[pidx] - dens_sorted[pidx]
                densrms = np.sqrt(np.mean(densanom ** 2))
                n2o = g * densrms / (Lto * np.mean(dens[pidx]))
            elif Nsq_method == "endpt":
                g = gsw.grav(lat, P[pidx].mean())
                ddens = dens_sorted[i1] - dens_sorted[i0]
                ddepth = depth[i1] - depth[i0]
                n2o = g * ddens / (ddepth * np.mean(dens[pidx]))
            else:
                raise ValueError("Nsq_method '{}' is not available.".format(Nsq_method))

            n2[pidx] = n2o

            # Flag beginning or end.
            if i0 == 0:
                ends_flag[pidx] = True
            if i1 == ndata - 1:
                ends_flag[pidx] = True

            # Flag small density or CT difference.
            if not use_ip:
                dq = q_sorted[i1] - q_sorted[i0]
                if dq < dnoise:
                    noise_flag[pidx] = True

            # Flag negative N squared.
            if n2o < 0:
                n2_flag[pidx] = True

            # Overturn ratio of Gargett & Garner
            Tdo = thorpe_disp[pidx]
            dzo = dz[pidx]
            L_tot = np.sum(dzo)
            L_neg = np.sum(dzo[Tdo < 0])
            L_pos = np.sum(dzo[Tdo > 0])
            Roo = np.minimum(L_neg / L_tot, L_pos / L_tot)
            Ro[pidx] = Roo
            # Don't raise flag if using intermediate profile because
            # the method leads to sorting weirdness that messes with
            # the overturn ratio.
            if Roo < Roc and not use_ip:
                Ro_flag[pidx] = True

        # Find and select data for this reference pressure range only.
        inbin = (P > pbins[idx_bin]) & (P < pbins[idx_bin + 1])

        # Fill flags.
        diag["noise_flag"][inbin] = noise_flag[inbin]
        diag["n2_flag"][inbin] = n2_flag[inbin]
        diag["ends_flag"][inbin] = ends_flag[inbin]
        diag["Ro_flag"][inbin] = Ro_flag[inbin]

        # Fill other diagnostics.
        diag["n2"][inbin] = n2[inbin]
        diag["Lt"][inbin] = Lt[inbin]
        diag["Ro"][inbin] = Ro[inbin]
        diag["thorpe_disp"][inbin] = thorpe_disp[inbin]
        diag["dens"][inbin] = dens[inbin]
        diag["dens_sorted"][inbin] = dens_sorted[inbin]
        diag["CT_sorted"][inbin] = CT_sorted[inbin]
        diag["SA_sorted"][inbin] = SA_sorted[inbin]

        if use_ip and not overturns_from_CT:
            diag["dens_ip"][inbin] = q_ip[inbin]

        if use_ip and overturns_from_CT:
            diag["CT_ip"][inbin] = q_ip[inbin]

    # Finally calculate epsilon for diagnostics, avoid nans, inf and negative n2.
    isgood = np.isfinite(diag["n2"]) & np.isfinite(diag["Lt"]) & ~diag["n2_flag"]
    diag["eps"][isgood] = (
        alpha ** 2 * diag["Lt"][isgood] ** 2 * diag["n2"][isgood] ** 1.5
    )

    # Use flags to get rid of bad overturns in basic output
    isbad = diag["noise_flag"] | diag["n2_flag"] | diag["Ro_flag"]
    eps = diag["eps"].copy()
    eps[isbad] = np.nan
    n2 = diag["n2"].copy()
    n2[isbad] = np.nan

    # Fill with background epsilon
    eps[np.isnan(eps)] = background_eps

    if return_diagnostics:
        return eps, n2, diag
    else:
        return eps, n2


def find_overturns(x, combine_gap=0):
    """Find the indices of unstable patches.

    Parameters
    ----------
    x : array_like 1D
        Profile of some quantity from which overturns can be detected
        e.g. temperature or density.
    combine_gap : float, optional
        Combine overturns that are separated by less than a given number of points.
        Default is 0.

    Returns
    -------
    idx_sorted : 1D numpy array
        Indices that sort the data x.
    idx_patches : (N, 2) numpy array
        Start and end indices of the overturns.

    """
    idx = np.arange(x.size, dtype=int)
    idx_sorted = np.argsort(x, kind="mergesort")
    idx_cumulative = np.cumsum(idx_sorted - idx)
    idx_patches = _consec_blocks(np.where(idx_cumulative > 0)[0], combine_gap)
    return idx_sorted, idx_patches


def intermediate_profile_topdown(x, acc, hinge):
    """Generate an intermediate profile starting at x[0] moving along the array.

    See Ferron et. al. 1998 and Gargett and Garner 2008.

    Parameters
    ----------
    x : array_like 1D
        Profile of some quantity that the intermediate profile method can be
        applied to e.g. temperature or density.
    acc : float, optional
        Accuracy parameter. The intermediate profile change in steps of acc.
    hinge : float, optional
        Intermediate profile values are equal to the hinge plus an integer
        multiple of acc. It should be kept constant across profiles.

    Returns
    -------
    xi : 1D numpy array
        Intermediate profile.

    """

    # Initialise.
    xi = np.zeros_like(x)
    n = np.fix((x[0] - hinge) / acc)
    xi[0] = hinge + n * acc

    # Step through profile.
    for i in range(len(x) - 1):
        n = np.fix((x[i + 1] - xi[i]) / acc)
        xi[i + 1] = xi[i] + n * acc

    return xi


def intermediate_profile(x, acc=5e-4, hinge=1000, kind="down"):
    """Generate an intermediate profile of some quantity.

    See Ferron et. al. 1998 and Gargett and Garner 2008.

    Parameters
    ----------
    x : array_like 1D
        Profile of some quantity that the intermediate profile method can be
        applied to e.g. temperature or density.
    acc : float, optional
        Accuracy parameter. The intermediate profile change in steps of acc.
    hinge : float, optional
        Intermediate profile values are equal to the hinge plus an integer multiple
        of acc. It should be kept constant across profiles.
    kind : string, optional
        Either 'up', 'down' or 'ave'. Default is ave. This argument determines
        whether the method is applied top down (x[0] to x[end]), bottom up
        (x[end] to [x[0]]) or the average of up and down.

    Returns
    -------
    xi : 1D numpy array
        Intermediate profile.

    """

    if not any(s in kind for s in ["up", "do", "av"]):
        raise ValueError("The 'kind' argument must be 'up', 'down' or 'ave'.")

    x = np.asarray(x)

    if "up" in kind:
        xf = np.flipud(x)
        xi = np.flipud(intermediate_profile_topdown(xf, acc, hinge))
    elif "do" in kind:
        xi = intermediate_profile_topdown(x, acc, hinge)
    elif "av" in kind:
        xf = np.flipud(x)
        xtd = intermediate_profile_topdown(x, acc, hinge)
        xbu = np.flipud(intermediate_profile_topdown(xf, acc, hinge))
        xi = (xtd + xbu) / 2.0

    return xi


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
