import numpy as np

from mixsea import overturn


# Use the ctd profile defined as fixture in conftest.py
def test_overturn(ctd_profile):
    eps, n2 = overturn.nan_eps_overturn(
        ctd_profile["depth"],
        ctd_profile["t"],
        ctd_profile["SP"],
        ctd_profile["lon"][0],
        ctd_profile["lat"][0],
    )
    assert np.nanmean(eps) < 1e-4


def test_overturn_const_s(ctd_profile):
    eps, n2 = overturn.nan_eps_overturn(
        ctd_profile["depth"],
        ctd_profile["t"],
        35,
        ctd_profile["lon"][0],
        ctd_profile["lat"][0],
    )
    assert np.nanmean(eps) < 1e-4


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

    assert np.all(overturn.contiguous_regions(condition) == idx)
