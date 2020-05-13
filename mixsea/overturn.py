import gsw
import numpy as np


def nan_eps_overturn(
    P, T, S, lon, lat, dnoise=5e-4, alpha_sq=0.9, background_eps=np.nan, use_ip=True, return_diagnostics=False,
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
        Salinity [PSU]
    lon : float
        Longitude of observation
    lat : float
        Latitude of observation
    dnoise : float, optional
        Noise level of density [kg/m^3]. Default is 5e-4.
    alpha_sq : float, optional
        Square of proportionality constant between Thorpe and Ozmidov scale.
        Default is 0.9.
    background_eps : float, optional
        Background epsilon where no overturn detected. Defaults to numpy.nan.
    use_ip : boolean, optional
        Sets whether to use the intermediate profile method. Default is True. If True, 
        the dnoise parameter is passed as the `accuracy' argument of the intermediate
        profile method. 
    return_diagnostics : dict, optional
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
    
    # Find non-NaNs
    notnan = np.isfinite(P) & np.isfinite(T) & np.isfinite(S)
    isnan = ~notnan
    
    if isnan.sum() == 0:  # If there are no NaNs then return. 
        return eps_overturn(P, T, S, lon, lat, dnoise, alpha_sq, background_eps, return_diagnostics, use_ip)
    
    eps = np.full_like(P, np.nan)
    n2 = np.full_like(P, np.nan)
    
    if return_diagnostics:
        eps[notnan], n2[notnan], diag = eps_overturn(P[notnan], T[notnan], S[notnan], lon, lat, dnoise, alpha_sq, background_eps, use_ip, return_diagnostics)
        
        # Replace nans in diagnostics if the size and shape seems right:
        Nnotnans = notnan.sum()
        for key in diag:
            if (np.size(diag[key]) == Nnotnans) & (np.ndim(diag[key]) == 1):
                ar = np.full_like(P, np.nan)
                ar[notnan] = diag[key]
                diag[key] = ar  # This will wipe out the old item.
        
        return eps, n2, diag
        
    else:
        eps[notnan], n2[notnan] = eps_overturn(P[notnan], T[notnan], S[notnan], lon, lat, dnoise, alpha_sq, background_eps, use_ip)
        return eps, n2
    


def eps_overturn(
    P, T, S, lon, lat, dnoise=5e-4, alpha_sq=0.9, background_eps=np.nan, use_ip=True, return_diagnostics=False,
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
    S : array-like
        Salinity [PSU]
    lon : float
        Longitude of observation
    lat : float
        Latitude of observation
    dnoise : float, optional
        Noise level of density [kg/m^3]. Default is 5e-4.
    alpha_sq : float, optional
        Square of proportionality constant between Thorpe and Ozmidov scale.
        Default is 0.9.
    background_eps : float, optional
        Background epsilon where no overturn detected. Defaults to numpy.nan.
    use_ip : boolean, optional
        Sets whether to use the intermediate profile method. Default is True. If True, 
        the dnoise parameter is passed as the `accuracy' argument of the intermediate
        profile method. 
    return_diagnostics : dict, optional
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
    # Estimate depth from pressure.
    depth = -gsw.z_from_p(P, lat)

    SA = gsw.SA_from_SP(S, T, lon, lat)
    CT = gsw.CT_from_t(SA, T, P)
    
    # Initialise diagnostics.
    diag = {}
    diag["eps"] = np.full_like(P, np.nan)
    diag["n2"] = np.full_like(P, np.nan)
    diag["Lt"] = np.full_like(P, np.nan)
    diag["dtdz"] = np.full_like(P, np.nan)
    diag["noise_flag"] = np.full_like(P, False, dtype=bool)

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
        # Calculate potential density
        dens = gsw.pot_rho_t_exact(SA, T, P, p_ref=refdi)

        if use_ip: # Create intermediate density profile
            dens_ip = intermediate_profile(dens, acc=dnoise, hinge=dens[0], kind="down")
            Is, patches = find_overturns(dens_ip)
        else:
            Is, patches = find_overturns(dens)

        # Calculate Thorpe displacements
        thorpe_displacement = depth[Is] - depth
        
        # Sort temperature and salinity based on the density sorting index
        # for calculating the buoyancy frequency
        SA_sorted = SA[Is]
        CT_sorted = CT[Is]
        dens_sorted = dens[Is]

        # Make sure there are any overturns.
        if np.any(patches):

            # Loop over detected overturns and calculate Thorpe Scales, N2
            # and dT/dz over the overturn region
            Lt = np.full_like(P, np.nan)
            N2 = np.full_like(P, np.nan)
            dTdz = np.full_like(P, np.nan)
            noise_flag = np.full_like(P, False, dtype=bool)

            for patch in patches:
                i0, i1 = patch[0], patch[1]
                idx = np.arange(i0, i1 + 1, 1)
                
                # If we are not using the intermediate profile method, need to check the density difference is greater than 
                # the noise threshold.
                if ~use_ip:
                    dens_diff = dens_sorted[i1+1] - dens_sorted[i0]
                    if dens_diff < dnoise:
                        noise_flag[idx] = True
                
                # Thorpe scale is the root mean square thorpe displacement. 
                Lt[idx] = np.sqrt(np.mean(np.square(thorpe_displacement[idx])))
                
                # Buoyancy frequency calculated over the overturn from sorted
                # profiles. Go beyond overturn.
                N2[idx], _ = gsw.Nsquared(
                    SA_sorted[[i0 - 1, i1 + 1]],
                    CT_sorted[[i0 - 1, i1 + 1]],
                    P[[i0 - 1, i1 + 1]],
                    lat,
                )

                # Fill depth range of the overturn with local temperature gradient
                if i0 > 0:
                    PTov = CT_sorted[i0 - 1 : i1 + 1]
                    zov = depth[i0 - 1 : i1 + 1]
                else:
                    PTov = CT_sorted[i0 : i1 + 1]
                    zov = depth[i0 : i1 + 1]
                local_dTdz = (np.min(PTov) - np.max(PTov)) / (np.max(zov) - np.min(zov))
                dTdz[idx] = local_dTdz

            # Find data for this reference depth range
            iz = np.squeeze(
                np.where(((P > refdi - dref / 2) & (P <= refdi + dref / 2)))
            )

            # Pick only data for this reference depth range
            diag["n2"][iz] = N2[iz]
            diag["Lt"][iz] = Lt[iz]
            diag["dtdz"][iz] = dTdz[iz]
            diag["noise_flag"][iz] = noise_flag[iz]
            
        # Finally calculate epsilon for diagnostics, avoid nans and inf.
        isfinite = np.isfinite(diag["n2"]) & np.isfinite(diag["Lt"])
        diag['eps'][isfinite] = alpha_sq * diag["Lt"][isfinite] ** 2 * diag["n2"][isfinite] ** 1.5
        
        # Use flags to get rid of bad overturns in basic output
        isbad = diag["noise_flag"]
        eps = diag['eps'].copy()
        eps[isbad] = np.nan
        n2 = diag['n2'].copy()
        n2[isbad] = np.nan
            
        # Fill with background epsilon
        eps[np.isnan(eps)] = background_eps

    if return_diagnostics:
        return eps, n2, diag
    else:
        return eps, n2
    

def find_overturns(x):
    """Find the indices of unstable patches.

    Parameters
    ----------
    x : array_like 1D
        Profile of some quantity from which overturns can be detected 
        e.g. temperature or density.

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
    idx_patches = _consec_blocks(np.where(idx_cumulative > 2)[0], combine_gap=1)
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
