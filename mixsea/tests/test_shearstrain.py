import numpy as np

from mixsea import shearstrain


# Use the ctd and ladcp profile defined as fixture in conftest.py
def test_shearstrain(ctd_profile, ladcp_profile):
    assert ctd_profile["lon"].shape == ctd_profile["lat"].shape
    # Center points of depth windows. Windows are half overlapping, i.e.
    # their size (300m) is double the spacing here (150m).
    window_size = 300
    dz = window_size / 2
    zbin = np.linspace(dz, dz * 40, num=40)
    # Wavenumber vector. Starts at wavenumber corresponding to window size.
    m = np.arange(
        start=2 * np.pi / window_size, stop=2 * np.pi / 10, step=2 * np.pi / window_size
    )
    # Wavenumber indices for integration. Shear is integrated from 300m to
    # 100m scales. Strain is integrated from 150m to 30m.
    m_include_sh = list(range(3))
    m_include_st = list(range(1, 10))

    (
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
    ) = shearstrain.shearstrain(
        ctd_profile["s"],
        ctd_profile["t"],
        ctd_profile["p"],
        ctd_profile["z"],
        ctd_profile["lat"],
        ctd_profile["lon"],
        ladcp_profile["uz"],
        ladcp_profile["vz"],
        ladcp_profile["z"],
        m=m,
        z_bin=zbin,
        m_include_sh=m_include_sh,
        m_include_st=m_include_st,
        ladcp_is_shear=True,
    )

    assert np.nanmean(np.log10(eps_shst)) < 0
