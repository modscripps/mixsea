import numpy as np


# Use the ctd profile defined as a fixture in conftest.py
def test_ctd(ctd_profile):
    assert np.mean(ctd_profile["lon"]) < 0


# Use the ladcp profile defined as a fixture in conftest.py
def test_ladcp(ladcp_profile):
    assert np.mean(ladcp_profile["depth"]) > 0
