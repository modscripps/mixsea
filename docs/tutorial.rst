=========
Tutorials
=========

For all tutorials, import mixsea and other libraries:

.. ipython:: python

    import mixsea as mx
    import numpy as np
    import matplotlib.pyplot as plt

Read example CTD/LADCP profiles
-------------------------------

Read the example CTD and LADCP profiles from the test suite using a little helper function:

.. ipython:: python

    ctd = mx.helpers.read_ctd_testfile()
    ladcp = mx.helpers.read_ladcp_testfile()

A quick overview plot of the data:

.. ipython:: python

    fig, ax = plt.subplots(nrows=1,
                           ncols=3,
                           figsize=(7, 3),
                           constrained_layout=True, 
                           sharey=True)
    ax[0].plot(ctd['t'], ctd['z']);
    ax[1].plot(ladcp['u'], ladcp['z']);
    ax[1].plot(ladcp['v'], ladcp['z']);
    ax[2].plot(ladcp['uz'], ladcp['z']);
    ax[0].invert_yaxis()
    @savefig ctd_test_profile.png width=7in 
    ax[0].set(ylabel='depth [m]');

This is cast 81 from the 2012 Samoan Passage cruise. The layer of cold
Antarctic Bottom Water flowing through the Samoan Passage shows up in
temperature, velocity and shear. See :cite:`voetetal15` for more information.


Thorpe scales / overturns
-------------------------



Shear/strain
------------

