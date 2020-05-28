import numpy as np

from mixsea import overturn


# Use the ctd profile defined as fixture in conftest.py
def test_overturn(ctd_profile):
    eps, n2 = overturn.nan_eps_overturn(
        ctd_profile["z"],
        ctd_profile["t"],
        ctd_profile["s"],
        ctd_profile["lon"][0],
        ctd_profile["lat"][0],
    )
    assert np.nanmean(eps) < 1e-4


def test_overturn_const_s(ctd_profile):
    eps, n2 = overturn.nan_eps_overturn(
        ctd_profile["z"],
        ctd_profile["t"],
        35,
        ctd_profile["lon"][0],
        ctd_profile["lat"][0],
    )
    assert np.nanmean(eps) < 1e-4
