import numpy as np
import pytest

from mixsea import shearstrain


# Use the ctd and ladcp profile defined as fixture in conftest.py
# Run the test both for adiabatic leveling and polynomial fit method.
# Use `nan_shearstrain` as the data contain NaNs.
@pytest.mark.parametrize("smooth", ["AL", "PF"])
def test_nan_shearstrain(ctd_profile, ladcp_profile, smooth):
    assert ctd_profile["lon"].shape == ctd_profile["lat"].shape
    # Center points of depth windows. Windows are half overlapping, i.e.
    # their size (300m) is double the spacing here (150m).
    window_size = 300
    dz = window_size / 2
    depth_bin = np.linspace(dz, dz * 40, num=40)
    # Wavenumber vector. Starts at wavenumber corresponding to window size.
    m = np.arange(
        start=2 * np.pi / window_size, stop=2 * np.pi / 10, step=2 * np.pi / window_size
    )
    # Wavenumber indices for integration. Shear is integrated from 300m to
    # 100m scales. Strain is integrated from 150m to 30m.
    m_include_sh = list(range(3))
    m_include_st = list(range(1, 10))

    (
        eps_shst,
        krho_shst,
        diag,
    ) = shearstrain.nan_shearstrain(
        ctd_profile["depth"],
        ctd_profile["t"],
        ctd_profile["SP"],
        ctd_profile["lon"],
        ctd_profile["lat"],
        ladcp_profile["uz"],
        ladcp_profile["vz"],
        ladcp_profile["depth"],
        m=m,
        depth_bin=depth_bin,
        window_size=window_size,
        m_include_sh=m_include_sh,
        m_include_st=m_include_st,
        ladcp_is_shear=True,
        smooth=smooth,
        return_diagnostics=True,
    )

    assert np.nanmean(np.log10(eps_shst)) < 0
    assert np.nanmean(np.log10(diag["eps_st"])) < 0


# Remove NaNs manually and use `shearstrain`.
def test_shearstrain(ctd_profile, ladcp_profile):
    ctd = ctd_profile
    notnan = nonan(ctd)
    depth = ctd["depth"][notnan]
    t = ctd["t"][notnan]
    SP = ctd["SP"][notnan]
    lon = ctd["lon"][0]
    lat = ctd["lat"][0]
    ladcp = ladcp_profile
    notnan = nonan(ladcp)
    uz = ladcp["uz"][notnan]
    vz = ladcp["vz"][notnan]
    zz = ladcp["depth"][notnan]

    # Center points of depth windows. Windows are half overlapping, i.e.
    # their size (300m) is double the spacing here (150m).
    window_size = 300
    dz = window_size / 2
    depth_bin = np.linspace(dz, dz * 40, num=40)
    # Wavenumber vector. Starts at wavenumber corresponding to window size.
    m = np.arange(
        start=2 * np.pi / window_size, stop=2 * np.pi / 10, step=2 * np.pi / window_size
    )
    # Wavenumber indices for integration. Shear is integrated from 300m to
    # 100m scales. Strain is integrated from 150m to 30m.
    m_include_sh = list(range(3))
    m_include_st = list(range(1, 10))

    (
        eps_shst,
        krho_shst,
        diag,
    ) = shearstrain.shearstrain(
        depth,
        t,
        SP,
        lon,
        lat,
        uz,
        vz,
        zz,
        m=m,
        depth_bin=depth_bin,
        window_size=window_size,
        m_include_sh=m_include_sh,
        m_include_st=m_include_st,
        ladcp_is_shear=True,
        smooth="PF",
        return_diagnostics=True,
    )

    assert np.nanmean(np.log10(eps_shst)) < 0
    assert np.nanmean(np.log10(diag["eps_st"])) < 0


def nonan(data):
    notnan = [np.isfinite(v) for k, v in data.items()]
    notnan = np.vstack(notnan)
    notnan = np.all(notnan, axis=0)
    return notnan
