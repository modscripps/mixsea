=======
History
=======

.. vx.x.x (unreleased)
.. ---------------------

.. Breaking changes
.. ~~~~~~~~~~~~~~~~

.. New Features
.. ~~~~~~~~~~~~

.. Bug fixes
.. ~~~~~~~~~

.. Documentation
.. ~~~~~~~~~~~~~

.. Internal Changes
.. ~~~~~~~~~~~~~~~~

v0.2.0 (unreleased)
---------------------

Breaking changes
~~~~~~~~~~~~~~~~
- Require the following minimum versions:
    - python>=3.9
    - numpy>=2
    - scipy>=1.6
    - gsw>=3.6

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
- Add zenodo doi to readme (:pull:`115`).
- Fix broken doi links (:pull:`117`).

Internal Changes
~~~~~~~~~~~~~~~~
- Modernize a number of internals:
    - Switch over to use `pyproject.toml` for python project definition and metadata.
    - Use `uv <https://docs.astral.sh/uv/>`_ as build backend.
    - Move code to `src/mixsea/`.
    - Move tests to root directory.
    - Remove pre-commit hooks for isort and black and add pre-commit hooks for `ruff <https://docs.astral.sh/ruff/>`_-check and ruff-format.
    - Remove test file with just a few examples (`tests/test_really_cool_feature.py`).
    - Update Github Workflow to work with uv and ruff.
    - Update Makefile to work with uv.
    - Remove `requirements.txt` and `requirements_docs.txt` as dependencies are now declared in `pyproject.toml`.

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
