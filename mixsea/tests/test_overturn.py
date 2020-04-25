import numpy as np

from mixsea import overturn


# Use the ctd profile defined as fixture in conftest.py
def test_overturn(ctd_profile):
    Lt, eps, k, n2, dtdz = overturn.eps_overturn(
        ctd_profile["p"],
        ctd_profile["z"],
        ctd_profile["t"],
        ctd_profile["s"],
        ctd_profile["lon"],
        ctd_profile["lat"],
    )
    assert np.nanmean(eps) < 1e-4
