import gsw
import numpy as np


def nan_eps_overturn(
    depth,
    t,
    SP,
    lon,
    lat,
    **kwargs,
):
    """
    Calculate turbulent dissipation based on the Thorpe scale method attempting to deal NaN values in the input data.
    It does this by removing all NaN values in the input profiles, then computes thorpe scales, then re-inserts NaNs
    at the end.

    See `eps_overturn` for more options.
    """
    depth = np.asarray(depth)
    t = np.asarray(t)
    SP = np.asarray(SP)

    # Find non-NaNs
    if SP.size == 1:
        SP = np.full_like(depth, SP)

    notnan = np.isfinite(depth) & np.isfinite(t) & np.isfinite(SP)

    isnan = ~notnan
    if isnan.sum() == 0:  # If there are no NaNs then return.
        return eps_overturn(depth, t, SP, lon, lat, **kwargs)

    eps = np.full_like(depth, np.nan)
    N2 = np.full_like(depth, np.nan)

    # Don't want to pass return_diagnostics twice.
    if "return_diagnostics" in kwargs:
        return_diagnostics = kwargs.pop("return_diagnostics")
    else:
        return_diagnostics = False

    eps[notnan], N2[notnan], diag = eps_overturn(
        depth[notnan],
        t[notnan],
        SP[notnan],
        lon,
        lat,
        return_diagnostics=True,
        **kwargs,
    )

    if return_diagnostics:
        # Replace nans in diagnostics if the size and shape seems right:
        Nnotnans = notnan.sum()
        for key in diag:
            if (np.size(diag[key]) == Nnotnans) & (np.ndim(diag[key]) == 1):
                ar = np.full_like(depth, np.nan)
                ar[notnan] = diag[key]
                diag[key] = ar  # This will wipe out the old item.

        return eps, N2, diag

    else:
        return eps, N2


def eps_overturn(
    depth,
    t,
    SP,
    lon,
    lat,
    dnoise=5e-4,
    alpha=0.95,
    Roc=0.2,
    background_eps=np.nan,
    use_ip=False,
    N2_method="teos",
    overturns_from_CT=False,
    return_diagnostics=False,
):
    """
    Calculate turbulent dissipation based on the Thorpe scale method. This function cannot handle
    NaNs in the input data, but there is another called `nan_eps_overturn' that attempts to.

    Parameters
    ----------
    depth : array-like
            Depth [m]
    t : array-like
            In-situ temperature [ITS90, °C]
    SP : float or array-like
            Salinity [PSU]. Can be a single constant value. This may be convenient if only temperature data
            are available.
    lon : float
            Longitude of observation
    lat : float
            Latitude of observation
    dnoise : float, optional
            Noise level of density [kg/m^3] or conservative temperature [°C], depending on overturns_from_CT. Default is 5e-4.
    alpha : float, optional
            Constant of proportionality between Thorpe and Ozmidov scale. Default is 0.95.
    Roc : float, optional
            Critical value of the overturn ratio Ro. An overturn will be considered
            noise if Ro < Roc.
    background_eps : float, optional
            Background epsilon where no overturn detected. Defaults to numpy.nan.
    use_ip : bool, optional
            Sets whether to use the intermediate profile method. Default is False. If True,
            the dnoise parameter is passed as the `accuracy' argument of the intermediate
            profile method.
    N2_method : string, optional
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
    eps : ndarray
            Turbulent dissipation [W/kg]
    N2 : ndarray
            Background stratification of each overturn detected [s^-2]
    diag : dict, optional
            Dictionary of diagnositc variables, set return with the `return_diagnostics' argument.


    """
    depth = np.asarray(depth)
    t = np.asarray(t)
    SP = np.asarray(SP)

    if not np.all(np.isclose(np.maximum.accumulate(depth), depth)):
        raise ValueError(
            "It appears that depth is not monotonically increasing, please fix."
        )

    if SP.size == 1:
        SP = np.full_like(depth, SP)

    if not (depth.size == t.size == SP.size):
        raise ValueError(
            "Input array sizes do not match. depth.size = {}, t.size = {}, SP.size = {}".format(
                depth.size, t.size, SP.size
            )
        )

    if not any(s == N2_method for s in ["teosp1", "teos", "endpt", "bulk"]):
        raise ValueError(
            "The 'N2_method' argument must be 'teosp1', 'teos', 'endpt' or 'bulk'."
        )

    ndata = depth.size

    # Estimate pressure from depth.
    p = gsw.p_from_z(-depth, lat)
    dz = 0.5 * (depth[2:] - depth[:-2])  # 'width' of each data point
    dz = np.hstack((dz[0], dz, dz[-1]))  # assume width of first and last data point

    SA = gsw.SA_from_SP(SP, t, lon, lat)
    CT = gsw.CT_from_t(SA, t, p)

    # Initialise arrays for diagnostic variables and flags.
    diag = {}
    diagvar = [
        "eps",
        "N2",
        "Lt",
        "thorpe_disp",
        "sidx",
        "dens",
        "dens_sorted",
        "Ro",
        "SA_sorted",
        "CT_sorted",
    ]
    for var in diagvar:
        diag[var] = np.full_like(depth, np.nan)

    if use_ip and not overturns_from_CT:
        diag["dens_ip"] = np.full_like(depth, np.nan)

    if use_ip and overturns_from_CT:
        diag["CT_ip"] = np.full_like(depth, np.nan)

    flagvar = ["noise_flag", "N2_flag", "ends_flag", "Ro_flag"]
    for var in flagvar:
        diag[var] = np.full_like(depth, False, dtype=bool)

    # Potential density is only meaningful near the reference pressure. For a deep profile
    # we may need to select several reference pressures. To do so, we find the pressure
    # bins that best contain the data
    pbinwidth = 1000.0  # In future we could have this as an argument.
    pbinmin = np.floor(p.min() / pbinwidth) * pbinwidth
    pbinmax = np.ceil(p.max() / pbinwidth) * pbinwidth
    pbins = np.arange(pbinmin, pbinmax + pbinwidth, pbinwidth)
    p_refs = 0.5 * (
        pbins[1:] + pbins[:-1]
    )  # Use mid point pressure as reference pressure.
    nbins = p_refs.size

    # Loop over pressure bins.
    for idx_bin in range(nbins):

        dens = gsw.pot_rho_t_exact(SA, t, p, p_ref=p_refs[idx_bin])

        if overturns_from_CT:
            # Temperature normally decreases towards the bottom which would mean the
            # find_overturns algorithm thinks the whole water column is unstable! Minus fixes that.
            q = -CT
        else:
            q = dens

        if use_ip:  # Create intermediate density profile
            q = intermediate_profile(
                q, acc=dnoise, hinge=1000, kind="down"
            )  # TODO: make hinge optional

        # --->> THORPE SCALES <<---
        (
            Lt,
            thorpe_disp,
            q_sorted,
            noise_flag,
            ends_flag,
            Ro,
            idx_patches,
            sidx,
        ) = thorpe_scale(depth, q, dnoise)

        # If there are no overturns, move on to the next pressure bin.
        if not np.any(idx_patches):
            continue

        # Thorpe displacements (by definition relative to initial locations, so unsorted)
        unsidx = np.argsort(sidx)
        thorpe_disp = (depth[sidx] - depth)[unsidx]

        # Sort other quantities based on the sorting indices.
        dens_sorted = dens[sidx]
        SA_sorted = SA[sidx]
        CT_sorted = CT[sidx]

        # Temporary arrays.
        N2 = np.full_like(depth, np.nan)
        N2_flag = np.full_like(depth, False, dtype=bool)
        Ro_flag = np.full_like(depth, False, dtype=bool)

        for patch in idx_patches:
            # Get patch indices.
            i0 = patch[0]
            i1 = patch[1]
            pidx = np.arange(i0, i1 + 1, 1)  # Need +1 for Python indexing

            Lto = np.unique(Lt[pidx])

            # Estimate the buoyancy frequency.
            if N2_method == "teos":
                N2o, _ = gsw.Nsquared(
                    SA_sorted[[i0, i1]],
                    CT_sorted[[i0, i1]],
                    p[[i0, i1]],
                    lat,
                )
            elif N2_method == "teosp1":
                # Go beyond overturn. Need to add 1 for this, unless end or beginning.
                addi = 0 if i1 == ndata - 1 else 1
                subi = 0 if i0 == 0 else 1

                N2o, _ = gsw.Nsquared(
                    SA_sorted[[i0 - subi, i1 + addi]],
                    CT_sorted[[i0 - subi, i1 + addi]],
                    p[[i0 - subi, i1 + addi]],
                    lat,
                )
            elif N2_method == "bulk":
                g = gsw.grav(lat, p[pidx].mean())
                densanom = dens[pidx] - dens_sorted[pidx]
                densrms = np.sqrt(np.mean(densanom ** 2))
                N2o = g * densrms / (Lto * np.mean(dens[pidx]))
            elif N2_method == "endpt":
                g = gsw.grav(lat, p[pidx].mean())
                ddens = dens_sorted[i1] - dens_sorted[i0]
                ddepth = depth[i1] - depth[i0]
                N2o = g * ddens / (ddepth * np.mean(dens[pidx]))
            else:  # May be redundent because of check at beginning of function.
                raise ValueError("N2_method '{}' is not available.".format(N2_method))

            N2[pidx] = N2o

            # Flag negative N squared.
            if N2o < 0:
                N2_flag[pidx] = True

            Roo = np.unique(Ro[pidx])

            if Roo < Roc:
                Ro_flag[pidx] = True

        # Find and select data for this reference pressure range only.
        inbin = (p > pbins[idx_bin]) & (p <= pbins[idx_bin + 1])

        # Fill flags.
        diag["noise_flag"][inbin] = noise_flag[inbin]
        diag["N2_flag"][inbin] = N2_flag[inbin]
        diag["ends_flag"][inbin] = ends_flag[inbin]
        diag["Ro_flag"][inbin] = Ro_flag[inbin]

        # Fill other diagnostics.
        diag["N2"][inbin] = N2[inbin]
        diag["Lt"][inbin] = Lt[inbin]
        diag["Ro"][inbin] = Ro[inbin]
        diag["thorpe_disp"][inbin] = thorpe_disp[inbin]
        diag["sidx"][inbin] = sidx[inbin]
        diag["dens"][inbin] = dens[inbin]
        diag["dens_sorted"][inbin] = dens_sorted[inbin]
        diag["CT_sorted"][inbin] = CT_sorted[inbin]
        diag["SA_sorted"][inbin] = SA_sorted[inbin]

        if use_ip and not overturns_from_CT:
            diag["dens_ip"][inbin] = q[inbin]

        if use_ip and overturns_from_CT:
            diag["CT_ip"][inbin] = q[inbin]

    # Finally calculate epsilon for diagnostics, avoid nans, inf and negative N2.
    isgood = np.isfinite(diag["N2"]) & np.isfinite(diag["Lt"]) & ~diag["N2_flag"]
    diag["eps"][isgood] = (
        alpha ** 2 * diag["Lt"][isgood] ** 2 * diag["N2"][isgood] ** 1.5
    )

    # Use flags to get rid of bad overturns in basic output
    isbad = diag["noise_flag"] | diag["N2_flag"] | diag["Ro_flag"]
    eps = diag["eps"].copy()
    eps[isbad] = np.nan
    N2 = diag["N2"].copy()
    N2[isbad] = np.nan

    # Fill with background epsilon
    eps[np.isnan(eps)] = background_eps

    if return_diagnostics:
        return eps, N2, diag
    else:
        return eps, N2


def thorpe_scale(depth, q, dnoise):
    """
    Estimate the Thorpe scale from unstable patches in a profile.

    Parameters
    ----------
    depth : array-like
            Depth [m] (negative if below sea surface)
    q : array-like
            Quantity from which Thorpe scales will be computed, e.g. density or temperature. If using
            temperature, consider multiplying by -1 to get around the fact that temperature generally
            decreases with depth.
    dnoise : float, optional
            Uncertainty or noise in q.

    Returns
    -------
    Lt : ndarray
            Thorpe scale [m]
    thorpe_disp : ndarray
            Thorpe displacement [m]
    q_sorted : ndarray
            q sorted to be monotonically increasing
    noise_flag : ndarray
            True if difference in q from top to bottom patch is less than dnoise
    ends_flag : ndarray
            True if a patch includes and end point
    Ro : ndarray
            Overturn ratio of Gargett & Garner.
    idx_patches : ndarray
            Indices of overturning patches, e.g. idx_patches[:, 0] are start indices and idx_patches[:, 1] are end indices.
    idx_sorted : ndarray
            Indices required to sort q so as to generate q_sorted.
    """

    depth = np.asarray(depth)
    q = np.asarray(q)

    if q[0] > q[-1]:
        raise ValueError("The entire profile is unstable, q[0] > q[-1].")

    if not np.all(np.isclose(np.maximum.accumulate(depth), depth)):
        raise ValueError(
            "It appears that depth is not monotonically increasing, please fix."
        )

    idx_sorted, idx_patches = find_overturns(q)

    ndata = depth.size

    # Thorpe displacements
    thorpe_disp = depth[idx_sorted] - depth

    q_sorted = q[idx_sorted]

    # Initialise arrays.
    Lt = np.full_like(depth, np.nan)
    Ro = np.full_like(depth, np.nan)
    noise_flag = np.full_like(depth, False, dtype=bool)
    ends_flag = np.full_like(depth, False, dtype=bool)

    dz = 0.5 * (depth[2:] - depth[:-2])  # 'width' of each data point
    dz = np.hstack((dz[0], dz, dz[-1]))  # assume width of first and last data point

    for patch in idx_patches:
        # Get patch indices.
        i0 = patch[0]
        i1 = patch[1]
        pidx = np.arange(i0, i1 + 1, 1)  # Need +1 for Python indexing

        # Thorpe scale is the root mean square thorpe displacement.
        Lto = np.sqrt(np.mean(np.square(thorpe_disp[pidx])))
        Lt[pidx] = Lto

        # Flag beginning or end.
        if i0 == 0:
            ends_flag[pidx] = True
        if i1 == ndata - 1:
            ends_flag[pidx] = True

        # Flag small difference.
        dq = q_sorted[i1] - q_sorted[i0]
        if dq < dnoise:
            noise_flag[pidx] = True

        # Overturn ratio of Gargett & Garner
        Tdo = thorpe_disp[pidx]
        dzo = dz[pidx]
        L_tot = np.sum(dzo)
        L_neg = np.sum(dzo[Tdo < 0])
        L_pos = np.sum(dzo[Tdo > 0])
        Roo = np.minimum(L_neg / L_tot, L_pos / L_tot)
        Ro[pidx] = Roo

    return Lt, thorpe_disp, q_sorted, noise_flag, ends_flag, Ro, idx_patches, idx_sorted


def find_overturns(q):
    """Find the indices of unstable patches by cumulatively summing the difference between
    sorted and unsorted indices of q.

    Parameters
    ----------
    q : array_like 1D
            Profile of some quantity from which overturns can be detected
            e.g. temperature or density.

    Returns
    -------
    idx_sorted : 1D ndarray
            Indices that sort the data q.
    idx_patches : (N, 2) ndarray
            Start and end indices of the overturns.

    """
    idx = np.arange(q.size, dtype=int)
    idx_sorted = np.argsort(q, kind="mergesort")
    idx_cumulative = np.cumsum(idx_sorted - idx)
    idx_patches = contiguous_regions(idx_cumulative > 0)
    return idx_sorted, idx_patches


def intermediate_profile_topdown(q, acc, hinge):
    """Generate an intermediate profile starting at q[0] moving along the array.

    See Ferron et. al. 1998 and Gargett and Garner 2008.

    Parameters
    ----------
    q : array_like 1D
            Profile of some quantity that the intermediate profile method can be
            applied to e.g. temperature or density.
    acc : float, optional
            Accuracy parameter. The intermediate profile change in steps of acc.
    hinge : float, optional
            Intermediate profile values are equal to the hinge plus an integer
            multiple of acc. It should be kept constant across profiles.

    Returns
    -------
    qi : 1D ndarray
            Intermediate profile.

    """

    # Initialise.
    qi = np.zeros_like(q)
    n = np.fix((q[0] - hinge) / acc)
    qi[0] = hinge + n * acc

    # Step through profile.
    for i in range(len(q) - 1):
        n = np.fix((q[i + 1] - qi[i]) / acc)
        qi[i + 1] = qi[i] + n * acc

    return qi


def intermediate_profile(q, acc=5e-4, hinge=1000, kind="down"):
    """Generate an intermediate profile of some quantity.

    See Ferron et. al. 1998 and Gargett and Garner 2008.

    Parameters
    ----------
    q : array_like 1D
            Profile of some quantity that the intermediate profile method can be
            applied to e.g. temperature or density.
    acc : float, optional
            Accuracy parameter. The intermediate profile change in steps of acc.
    hinge : float, optional
            Intermediate profile values are equal to the hinge plus an integer multiple
            of acc. It should be kept constant across profiles.
    kind : string, optional
            Either 'up', 'down' or 'ave'. Default is ave. This argument determines
            whether the method is applied top down (q[0] to q[end]), bottom up
            (q[end] to [q[0]]) or the average of up and down.

    Returns
    -------
    qi : 1D ndarray
            Intermediate profile.

    """

    if not any(s in kind for s in ["up", "do", "av"]):
        raise ValueError("The 'kind' argument must be 'up', 'down' or 'ave'.")

    q = np.asarray(q)

    if "up" in kind:
        qf = np.flipud(q)
        qi = np.flipud(intermediate_profile_topdown(qf, acc, hinge))
    elif "do" in kind:
        qi = intermediate_profile_topdown(q, acc, hinge)
    elif "av" in kind:
        qf = np.flipud(q)
        qtd = intermediate_profile_topdown(q, acc, hinge)
        qbu = np.flipud(intermediate_profile_topdown(qf, acc, hinge))
        qi = (qtd + qbu) / 2.0

    return qi


def contiguous_regions(condition):
    """Finds the indices of contiguous True regions in a boolean array.

    Parameters
    ----------
    condition : array_like
            Array of boolean values.

    Returns
    -------
    idx : ndarray
            Array of indices demarking the start and end of contiguous True regions in condition.
            Shape is (N, 2) where N is the number of regions.

    Notes
    -----
    Modified from stack overflow: https://stackoverflow.com/a/4495197

    """

    d = np.diff(condition)
    (idx,) = d.nonzero()

    # We need to start things after the change in "condition". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size]  # Edit

    # Reshape the result into two columns
    idx.shape = (-1, 2)
    return idx
