import numpy as np

from mixsea import overturn


# Use the ctd profile defined as fixture in conftest.py
def test_nan_eps_overturn(ctd_profile):
    # This also tests eps_overturn

    # First check that we get something out with default parameters
    eps, N2 = overturn.nan_eps_overturn(
        ctd_profile["depth"],
        ctd_profile["t"],
        ctd_profile["SP"],
        lon=ctd_profile["lon"][0],
        lat=ctd_profile["lat"][0],
    )
    assert np.nanmean(eps) < 1e-4
    assert (np.isfinite(eps) == np.isfinite(N2)).all()
    assert (N2[np.isfinite(N2)] < 1e-1).all()

    # Perform same checks as above for constant salinity input
    eps, N2 = overturn.nan_eps_overturn(
        ctd_profile["depth"],
        ctd_profile["t"],
        35,
        lon=ctd_profile["lon"][0],
        lat=ctd_profile["lat"][0],
    )
    assert np.nanmean(eps) < 1e-4
    assert (np.isfinite(eps) == np.isfinite(N2)).all()
    assert (N2[np.isfinite(N2)] < 1e-1).all()


def test_eps_overturn():
    # Generate synthetic profile with a single overturn
    depth = np.linspace(1, 50, 10000)
    SP = np.linspace(26, 29, 10000)
    t = np.linspace(4, 5, 10000)
    lon = -60
    lat = 60
    h = 10  # overturn height

    # Create overturn
    i0, i1 = np.searchsorted(depth, [20, 20 + h])
    SPo = SP.copy()
    to = t.copy()
    SPo[i0:i1] = SP[i1 - 1 : i0 - 1 : -1]
    to[i0:i1] = t[i1 - 1 : i0 - 1 : -1]

    # Estimate the thermodynamics, as it should be done in eps_overturn
    p = overturn.gsw.p_from_z(-depth, lat)
    SA = overturn.gsw.SA_from_SP(SP, p, lon, lat)
    CT = overturn.gsw.CT_from_t(SA, t, p)
    dens = overturn.gsw.pot_rho_t_exact(SA, t, p, 500.0)
    N2p, _ = overturn.gsw.Nsquared(SA, CT, p, lat)  # Almost constant...

    SAo = overturn.gsw.SA_from_SP(SPo, p, lon, lat)
    denso = overturn.gsw.pot_rho_t_exact(SAo, to, p, 500.0)

    eps, N2, diag = overturn.eps_overturn(
        depth, to, SPo, lon, lat, alpha=1.0, return_diagnostics=True
    )

    inoverturn = eps > 0

    Lt_true = h / np.sqrt(3)
    N2_true = N2p.mean()
    eps_true = Lt_true**2 * N2_true**1.5

    assert np.isclose(diag["dens"], denso).all()
    assert np.isclose(np.unique(diag["N2"][inoverturn])[0], N2_true, atol=1e-5)
    assert np.isclose(np.unique(diag["Lt"][inoverturn])[0], Lt_true, atol=5e-2)
    assert np.isclose(np.unique(diag["eps"][inoverturn])[0], eps_true, atol=1e-5)

    # Now check bulk method
    eps, N2, diag = overturn.eps_overturn(
        depth, to, SPo, lon, lat, alpha=1.0, N2_method="bulk", return_diagnostics=True
    )
    inoverturn = eps > 0

    assert np.isclose(np.unique(diag["N2"][inoverturn])[0], N2_true, atol=1e-5)
    assert np.isclose(np.unique(diag["Lt"][inoverturn])[0], Lt_true, atol=5e-2)
    assert np.isclose(np.unique(diag["eps"][inoverturn])[0], eps_true, atol=1e-5)

    # Now check endpt method
    eps, N2, diag = overturn.eps_overturn(
        depth, to, SPo, lon, lat, alpha=1.0, N2_method="endpt", return_diagnostics=True
    )
    inoverturn = eps > 0

    assert np.isclose(np.unique(diag["N2"][inoverturn])[0], N2_true, atol=1e-5)
    assert np.isclose(np.unique(diag["Lt"][inoverturn])[0], Lt_true, atol=5e-2)
    assert np.isclose(np.unique(diag["eps"][inoverturn])[0], eps_true, atol=1e-5)

    # Now check linear equation of state
    eps, N2, diag = overturn.eps_overturn(
        depth,
        to,
        SPo,
        lon,
        lat,
        alpha=1.0,
        N2_method="endpt",
        EOS="linear",
        return_diagnostics=True,
    )
    inoverturn = eps > 0
    N2_true = -9.82 * (
        -2e-4 * np.gradient(t, -depth).mean() + 7e-4 * np.gradient(SP, -depth).mean()
    )
    eps_true = Lt_true**2 * N2_true**1.5

    assert np.isclose(np.unique(diag["N2"][inoverturn])[0], N2_true, atol=1e-5)
    assert np.isclose(np.unique(diag["Lt"][inoverturn])[0], Lt_true, atol=5e-2)
    assert np.isclose(np.unique(diag["eps"][inoverturn])[0], eps_true, atol=1e-5)

    eps, N2, diag = overturn.eps_overturn(
        depth,
        to,
        SPo,
        lon,
        lat,
        alpha=1.0,
        N2_method="bulk",
        EOS="linear",
        return_diagnostics=True,
    )
    inoverturn = eps > 0

    assert np.isclose(np.unique(diag["N2"][inoverturn])[0], N2_true, atol=1e-5)
    assert np.isclose(np.unique(diag["Lt"][inoverturn])[0], Lt_true, atol=5e-2)
    assert np.isclose(np.unique(diag["eps"][inoverturn])[0], eps_true, atol=1e-5)


def test_thorpe_scale():
    dz = 0.1
    z0 = -2200.0
    H = 100.0
    h = 50.0
    z = np.arange(z0 + H - dz / 2.0, z0 - H - dz / 2.0, -dz)
    rho_0 = 1041.4
    delta = 1e-4

    top = ((z0 + H) >= z) & (z >= (z0 + h))
    middle = ((z0 + h) > z) & (z > (z0 - h))
    bottom = ((z0 - h) >= z) & (z >= (z0 - H))

    # Thorpe 77 profile
    rho = np.zeros_like(z)
    rho[top] = rho_0 * (1.0 - delta)
    rho[middle] = rho_0 * (
        1.0 + delta * np.sin(3.0 * np.pi * (z[middle] - z0) / (2.0 * h))
    )
    rho[bottom] = rho_0 * (1.0 + delta)
    rho -= (z - z0) * 1.0e-13  # gets rid of spurious displacements in isothermal layers

    # Thorpe 77 analytical solution
    rho_s = np.zeros_like(z)
    rho_s[top] = rho[top]
    rho_s[middle] = rho_0 * (1.0 - delta * np.sin(np.pi * (z[middle] - z0) / (2.0 * h)))
    rho_s[bottom] = rho[bottom]

    Lt_analytical = 4.0 * h / (3.0 * np.sqrt(3.0))

    (
        Lt,
        thorpe_disp,
        q_sorted,
        noise_flag,
        ends_flag,
        Ro,
        patches,
        idx_sorted,
    ) = overturn.thorpe_scale(-z, rho, 0)

    inoverturn = thorpe_disp != 0.0

    assert np.isclose(q_sorted, rho_s).all()
    assert np.isclose(np.unique(Lt[inoverturn])[0], Lt_analytical)
    assert np.isclose(np.unique(Ro[inoverturn])[0], 0.5, atol=0.02)


def test_contiguous_regions():
    condition = np.array(
        [True, True, False, False, True, True, True, False, True, True]
    )
    idx = np.array([[0, 2], [4, 7], [8, 10]])

    assert (overturn.contiguous_regions(condition) == idx).all()


def test_pot_rho_linear():
    rho0 = 1025.0
    # With these inputs, the output should be rho0.
    pot_rho = overturn.pot_rho_linear(35, 15, rho0=1025, a=2e-4, b=7e-4, SP0=35, t0=15)
    assert np.isclose(rho0, pot_rho)


def test_find_overturns():
    q = np.array([1, 2, 3, 4, 6, 5, 7, 8, 9, 12, 11, 10, 13])
    idx_sorted, idx_patches = overturn.find_overturns(q)
    assert (idx_sorted == np.array([0, 1, 2, 3, 5, 4, 6, 7, 8, 11, 10, 9, 12])).all()
    assert (idx_patches == np.array([[4, 5], [9, 11]])).all()


def test_intermediate_profile_topdown():
    q = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 2.75, 3.0, 3.001])
    acc = 1.0
    hinge = 0.0
    qi = overturn.intermediate_profile_topdown(q, acc, hinge)
    assert (qi == np.array([0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0])).all()


def test_intermediate_profile():
    q = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 2.75, 3.0, 3.001])
    acc = 1.0
    hinge = 0.0
    kind = "down"
    qi = overturn.intermediate_profile(q, acc, hinge, kind)
    assert (qi == np.array([0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0])).all()
    kind = "up"
    qi = overturn.intermediate_profile(q, acc, hinge, kind)
    assert (qi == np.array([0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0])).all()
    kind = "av"
    qi = overturn.intermediate_profile(q, acc, hinge, kind)
    assert (qi == np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 2.5, 3.0, 3.0])).all()
