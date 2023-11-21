.. currentmodule:: mixsea

#############
API reference
#############

This page provides an auto-generated summary of mixsea's API. For more details
and examples, refer to the relevant chapters in the main part of the
documentation.

Overturn-based parameterization
===============================

.. autosummary::
    :toctree: generated/

    overturn.eps_overturn
    overturn.nan_eps_overturn
    overturn.thorpe_scale

Shear/strain-based parameterization
===================================

.. autosummary::
    :toctree: generated/

    shearstrain.shearstrain
    shearstrain.nan_shearstrain
    shearstrain.wavenumber_vector
    shearstrain.latitude_correction
    shearstrain.strain_polynomial_fits
    shearstrain.strain_adiabatic_leveling
    shearstrain.gm_shear_variance
    shearstrain.gm_strain_variance
    shearstrain.find_cutoff_wavenumber
    shearstrain.diffusivity
    shearstrain.aspect_ratio_correction_shst
    shearstrain.aspect_ratio_correction_st

Stratification
==============

.. autosummary::
    :toctree: generated/

    nsq.adiabatic_leveling

Other functions
===============

.. autosummary::
    :toctree: generated/

    helpers.psd
    helpers.read_ctd_testfile
    helpers.read_ladcp_testfile 
