=======
History
=======

v0.1.3 (unreleased)
---------------------

.. Breaking changes
.. ~~~~~~~~~~~~~~~~
    
New Features
~~~~~~~~~~~~
- Add shear-only solution to `diag` output in finescale parameterization (:pull:`119`).
  By `Gunnar Voet <https://github.com/gunnarvoet>`_.

.. Bug fixes
.. ~~~~~~~~~

Documentation
~~~~~~~~~~~~~
- Fix an error in stairs plot in shearstrain notebook (:pull:`112`). 
  :issue:`111` by `Ole Pinner <https://github.com/opinner>`_.
- Run `black` formatting on all jupyter notebooks (:pull:`113`). 
  By `Gunnar Voet <https://github.com/gunnarvoet>`_.

.. Internal Changes
.. ~~~~~~~~~~~~~~~~


v0.1.2 (2023-11-21)
---------------------

Breaking changes
~~~~~~~~~~~~~~~~
- In overturn.eps_overturn the argument overturns_from_CT was renamed to overturns_from_t (:pull:`97`). 
  By `Jesse Cusack <https://github.com/jessecusack>`_.
    
New Features
~~~~~~~~~~~~
- Linear equation of state option added to the overturn module (:pull:`97`) as well as a few other tweaks to the eps_overturn function, including:

    - making latitude and longitude arguments optional
    - providing an argument for the pressure bin width used in the potential density calculation
    - removing unnecessary and/or meaningless diagnostics 
  By `Jesse Cusack <https://github.com/jessecusack>`_.

Bug fixes
~~~~~~~~~
- Fix frequency shift bug in psd (:pull:`105`). 
  By `Gunnar Voet <https://github.com/gunnarvoet>`_.

Documentation
~~~~~~~~~~~~~
- Explanation of the linear equation of state (:pull:`97`).
  By `Jesse Cusack <https://github.com/jessecusack>`_.


Internal Changes
~~~~~~~~~~~~~~~~
- Many unit tests for the overturn module were added (:pull:`97`).
  By `Jesse Cusack <https://github.com/jessecusack>`_.


v0.1.1 (2022-05-12)
---------------------

This release brings lots of additions to the documentation and some other minor additions. We haven't gotten into the routine of adding changes to this file, so the notes below do not reflect all changes of this release.

Bug fixes
~~~~~~~~~
- Fix an indexing bug in the shear/strain parameterization (:pull:`80`).
  By `Henri Drake <https://github.com/hdrake>`_.


Documentation
~~~~~~~~~~~~~
- Lots of additions to the documentation (:pull:`80`).
  By `Henri Drake <https://github.com/hdrake>`_ and `Jesse Cusack <https://github.com/jessecusack>`_.


Internal Changes
~~~~~~~~~~~~~~~~
- Simplify overturn helper functions and unify variable names. (:pull:`93`).
  By `Jesse Cusack <https://github.com/jessecusack>`_.
- Add Thorpe scale function with unit test (:pull:`79`).
  By `Henri Drake <https://github.com/hdrake>`_ and `Jesse Cusack <https://github.com/jessecusack>`_.


0.1.0 (2020-06-02)
------------------

* Initial release.
