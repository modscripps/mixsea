.. .. highlight:: shell

************
Contributing
************

.. contents:: Table of contents:
    :local:

.. note::

  This document is based on ideas and guidelines from the
  `xarray Contributing Guide 
  <http://xarray.pydata.org/en/stable/contributing.html>`_,
  which in turn builds largely upon the
  `Pandas Contributing Guide 
  <http://pandas.pydata.org/pandas-docs/stable/contributing.html>`_ .
  Some of the guidelines are more aspirational than practical at this point and constitute an attempt to learn how to build a proper open source software repository.

Where to start?
===============

All contributions, bug reports, bug fixes, documentation improvements,
enhancements, and ideas are welcome.

If you are brand new to *mixsea* or open-source development, we recommend going
through the `GitHub "issues" tab <https://github.com/modscripps/mixsea/issues>`_
to find issues that interest you.
There are a number of issues listed under
`Documentation <https://github.com/modscripps/mixsea/issues?q=is%3Aissue+is%3Aopen+label%3Adocumentation>`_ and 
`good first issue <https://github.com/modscripps/mixsea/issues?q=is%3Aopen+is%3Aissue+label%3A%22good+first+issue%22>`_
where you could start out. Once you've found an interesting issue, you can return
here to get your development environment setup.

.. _contributing.bug_reports:

Bug reports and enhancement requests
====================================

Bug reports are an important part of making *mixsea* more stable. Having a complete bug
report will allow others to reproduce the bug and provide insight into fixing. See
`this stackoverflow article <https://stackoverflow.com/help/mcve>`_ for tips on
writing a good bug report.

Trying the bug-producing code out on the *main* branch is often a worthwhile exercise
to confirm the bug still exists. It is also worth searching existing bug reports and
pull requests to see if the issue has already been reported and/or fixed.

Bug reports must:

.. TODO: need to include a defined test data file for users in these examples.

#. Include a short, self-contained Python snippet reproducing the problem.
   You can format the code nicely by using `GitHub Flavored Markdown
   <http://github.github.com/github-flavored-markdown/>`_::

      ```python
      >>> import mixsea
      ...
      ```

.. 
  #. Include the full version string of *mixsea* and its dependencies. You can use the built in function::

      >>> import mixsea as ctd
      >>> ctd.show_versions()

#. Explain why the current behavior is wrong/not desired and what you expect instead.

The issue will then show up to the *mixsea* community and be open to comments/ideas
from others.

.. _contributing.github:

Working with the code
=====================

Now that you have an issue you want to fix, enhancement to add, or documentation
to improve, you need to learn how to work with GitHub and the *mixsea* code base.

.. _contributing.version_control:

Version control, Git, and GitHub
--------------------------------

To the new user, working with Git is one of the more daunting aspects of contributing
to *mixsea*.  It can very quickly become overwhelming, but sticking to the guidelines
below will help keep the process straightforward and mostly trouble free.  As always,
if you are having difficulties please feel free to ask for help.

The code is hosted on `GitHub <https://www.github.com/modscripps/mixsea>`_. To contribute you will need to sign up for a `free GitHub account
<https://github.com/signup/free>`_. We use `Git <http://git-scm.com/>`_ for
version control to allow many people to work together on the project.

Some great resources for learning Git:

* the `GitHub help pages <http://help.github.com/>`_.
* the `NumPy's documentation <http://docs.scipy.org/doc/numpy/dev/index.html>`_.
* Matthew Brett's `Pydagogue <http://matthew-brett.github.com/pydagogue/>`_.

Getting started with Git
------------------------

`GitHub has instructions <http://help.github.com/set-up-git-redirect>`__ for installing git,
setting up your SSH key, and configuring git.  All these steps need to be completed before
you can work seamlessly between your local repository and GitHub.

.. _contributing.forking:

Forking
-------

You will need your own fork to work on the code. Go to the `mixsea project
page <https://github.com/modscripps/mixsea>`_ and hit the ``Fork`` button. You will
want to clone your fork to your machine::

    git clone https://github.com/your-user-name/mixsea.git
    cd mixsea
    git remote add upstream https://github.com/modscripps/mixsea.git

This creates the directory `mixsea` and connects your repository to
the upstream (main project) *mixsea* repository.

.. _contributing.dev_env:

Creating a development environment
----------------------------------

To test out code changes, you'll need to build *mixsea* from source, which
requires a Python environment.

.. _contributiong.dev_python:

Creating a Python Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before starting any development, you'll need to create an isolated mixsea
development environment:

- Install either `Anaconda <https://www.anaconda.com/download/>`_ or `miniconda
  <https://conda.io/miniconda.html>`_
- Make sure your conda is up to date (``conda update conda``)
- Make sure that you have :ref:`cloned the repository <contributing.forking>`
- ``cd`` to the *mixsea* source directory

We'll now kick off a two-step process:

1. Install the build dependencies
2. Build and install mixsea

.. code-block:: none

   # Create and activate the build environment
   # This is for Linux and MacOS. On Windows, use py37-windows.yml instead.
   conda env create -f environment.yml

   conda activate mixsea

   # or with older versions of Anaconda:
   source activate mixsea

   # Build and install mixsea
   pip install -e .

At this point you should be able to import *mixsea* from your locally built version::

   $ python  # start an interpreter
   >>> import mixsea
   >>> mixsea.__version__
   '0.1.0'

This will create the new environment, and not touch any of your existing environments,
nor any existing Python installation.

To view your environments::

      conda info -e

To return to your root environment::

      conda deactivate

See the full conda docs `here <http://conda.pydata.org/docs>`__.

Creating a branch
-----------------

You want your main branch to reflect only production-ready code, so create a
feature branch for making your changes. For example::

    git branch shiny-new-feature
    git checkout shiny-new-feature

The above can be simplified to::

    git checkout -b shiny-new-feature

This changes your working directory to the shiny-new-feature branch.  Keep any
changes in this branch specific to one bug or feature so it is clear
what the branch brings to *mixsea*. You can have many "shiny-new-features"
and switch in between them using the ``git checkout`` command.

To update this branch, you need to retrieve the changes from the main branch::

    git fetch upstream
    git rebase upstream/main

This will replay your commits on top of the latest *mixsea* git main.  If this
leads to merge conflicts, you must resolve these before submitting your pull
request.  If you have uncommitted changes, you will need to ``git stash`` them
prior to updating.  This will effectively store your changes and they can be
reapplied after updating.

.. _contributing.documentation:

Contributing to the documentation
=================================

If you're not the developer type, contributing to the documentation is still of
huge value. The documentation is written in **reStructuredText**, which is almost like writing
in plain English, and built using `Sphinx <http://sphinx-doc.org/>`__. The
Sphinx Documentation has an excellent `introduction to reST
<http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`__.
Review the Sphinx docs to perform more complex changes to the documentation as well.

Some other important things to know about the docs:

- The *mixsea* documentation consists of two parts: the docstrings in the code
  itself and the docs in this folder ``mixsea/docs/``.

  The docstrings are meant to provide a clear explanation of the usage of the
  individual functions, while the documentation in this folder consists of
  tutorial-like overviews per topic together with some other information
  (what's new, installation, etc).

- The docstrings follow the **Numpy Docstring Standard**, which is used widely
  in the Scientific Python community. This standard specifies the format of
  the different sections of the docstring. See `this document
  <https://github.com/numpy/numpy/blob/main/doc/HOWTO_DOCUMENT.rst.txt>`_
  for a detailed explanation, or look at some of the existing functions to
  extend it in a similar manner.

- The tutorials make heavy use of the `ipython directive
  <http://matplotlib.org/sampledoc/ipython_directive.html>`_ sphinx extension.
  This directive lets you put code in the documentation which will be run
  during the doc build. For example::

      .. ipython:: python

          x = 2
          x**3

  will be rendered as::

      In [1]: x = 2

      In [2]: x**3
      Out[2]: 8

  Almost all code examples in the docs are run (and the output saved) during the
  doc build. This approach means that code examples will always be up to date,
  but it does make the doc building a bit more complex

- The documentation includes jupyer notebooks via the `nbsphinx <https://nbsphinx.readthedocs.io>`__ exension. Make sure to clean all notebook output before adding and committing changes, this removes any images and other output and makes diffs much more easy to read.

- Our API documentation in ``doc/api.rst`` houses the auto-generated documentation from the docstrings. Every method should be included in ``api.rst``, else Sphinx will emit a warning.


Building the documentation
==========================

Navigate to your root ``mixsea/`` directory in the console and run::

    make docs

Then you can find the HTML output in the folder ``mixsea/docs/_build/html/``.
The newly built site will also automatically open in your browser.

Alternatively, you can also run::

    make servedocs

This will watch for changes in the source files of the documentation and
rebuild whenever it detects a change.

.. _contributing.code:

Contributing to the code base
=============================

.. contents:: Code Base:
   :local:


Code standards
--------------

Writing good code is not just about what you write. It is also about *how* you
write it. During :ref:`Continuous Integration <contributing.ci>` testing, several
tools will be run to check your code for stylistic errors.
Generating any warnings will cause the test to fail.
Thus, good style is a requirement for submitting code to *mixsea*.

In addition, because other people may use our library, it is important that we
do not make sudden changes to the code that could have the potential to break
a lot of user code as a result, that is, we need it to be as *backwards compatible*
as possible to avoid mass breakages.

.. _code.formatting:

Code Formatting
~~~~~~~~~~~~~~~

mixsea uses several tools to ensure a consistent code format throughout the project:

- `Black <https://black.readthedocs.io/en/stable/>`_ for standardized code formatting
- `Flake8 <http://flake8.pycqa.org/en/latest/>`_ for general code quality
- `isort <https://github.com/timothycrosley/isort>`_ for standardized order in imports.
  See also `flake8-isort <https://github.com/gforcada/flake8-isort>`_.

``pip``::

   pip install black flake8 isort mypy

and then run from the root of the `mixsea` repository::

   isort -rc .
   black -t py36 .
   flake8

to auto-format your code. Additionally, many editors have plugins that will
apply ``black`` as you edit files.

Optionally, you may wish to setup `pre-commit hooks <https://pre-commit.com/>`_
to automatically run all the above tools every time you make a git commit. This
can be done by installing ``pre-commit``::

   pip install pre-commit

and then running::

   pre-commit install

from the root of the mixsea repository. You can skip the pre-commit checks
with ``git commit --no-verify``.


.. _contributing.ci:

Testing With Continuous Integration
-----------------------------------

The *mixsea* test suite runs automatically the
`Travis CI <https://travis-ci.com/github/modscripps/mixsea>`__,
continuous integration service, once your pull request is submitted.

A pull-request will be considered for merging when you have an all 'green' build. If any tests are failing, then you will get a red 'X', where you can click through to see the individual failed tests.

.. note::

   Each time you push to your PR branch, a new run of the tests will be
   triggered on the CI.


.. _contributing.tdd:

Test-driven development/code writing
------------------------------------

*mixsea* is serious about testing and strongly encourages contributors to embrace
`test-driven development (TDD) <http://en.wikipedia.org/wiki/Test-driven_development>`_.
This development process "relies on the repetition of a very short development cycle:
first the developer writes an (initially failing) automated test case that defines a desired
improvement or new function, then produces the minimum amount of code to pass that test."
So, before actually writing any code, you should write your tests.  Often the test can be
taken from the original GitHub issue.  However, it is always worth considering additional
use cases and writing corresponding tests.

.. Adding tests is one of the most common requests after code is pushed to *mixsea*.  Therefore, it is worth getting in the habit of writing tests ahead of time so this is never an issue.

Like many packages, *mixsea* uses `pytest
<http://doc.pytest.org/en/latest/>`_ and the convenient
extensions in `numpy.testing
<http://docs.scipy.org/doc/numpy/reference/routines.testing.html>`_.

Writing tests
~~~~~~~~~~~~~

All tests should go into the ``tests`` subdirectory of the specific package.
This folder contains current examples of tests, and we suggest looking to these for
inspiration.

Here is an example of a self-contained set of tests that illustrate multiple
features that we like to use.

- functional style: tests are like ``test_*`` and *only* take arguments that are either
  fixtures or parameters
- ``pytest.mark`` can be used to set metadata on test functions, e.g. ``skip`` or ``xfail``.
- using ``parametrize``: allow testing of multiple cases
- to set a mark on a parameter, ``pytest.param(..., marks=...)`` syntax should be used
- ``fixture``, code for object construction, on a per-test basis
- using bare ``assert`` for scalars and truth-testing


We would name this file ``test_really_cool_feature.py`` and put in an appropriate place in the
``mixsea/tests/`` structure.

.. code-block:: python

    import pytest
    import numpy as np


    @pytest.mark.parametrize('dtype', ['int8', 'int16', 'int32', 'int64'])
    def test_dtypes(dtype):
        assert str(np.dtype(dtype)) == dtype


    @pytest.mark.parametrize('dtype', ['float32',
                             pytest.param('int16', marks=pytest.mark.skip),
                             pytest.param('int32', marks=pytest.mark.xfail(
                                reason='to show how it works'))])
    def test_mark(dtype):
        assert str(np.dtype(dtype)) == 'float32'


    @pytest.fixture
    def fake_data():
        return np.array([1, 2, 3])


    @pytest.fixture(params=['int8', 'int16', 'int32', 'int64'])
    def dtype(request):
        return request.param


    def test_series(fake_data, dtype):
        result = fake_data.astype(dtype)
        assert result.dtype == dtype



A test run of this yields

.. code-block:: shell

    (mixsea) $ pytest mixsea/tests/test_really_cool_feature.py -v
    ========================================= test session starts =====================
    platform darwin -- Python 3.8.2, pytest-5.4.1, py-1.8.1, pluggy-0.13.1 --
    cachedir: .pytest_cache
    plugins: datadir-1.3.1
    collected 11 items

    mixsea/tests/test_really_cool_feature.py::test_dtypes[int8] PASSED          [  9%]
    mixsea/tests/test_really_cool_feature.py::test_dtypes[int16] PASSED         [ 18%]
    mixsea/tests/test_really_cool_feature.py::test_dtypes[int32] PASSED         [ 27%]
    mixsea/tests/test_really_cool_feature.py::test_dtypes[int64] PASSED         [ 36%]
    mixsea/tests/test_really_cool_feature.py::test_mark[float32] PASSED         [ 45%]
    mixsea/tests/test_really_cool_feature.py::test_mark[int16] SKIPPED          [ 54%]
    mixsea/tests/test_really_cool_feature.py::test_mark[int32] XFAIL            [ 63%]
    mixsea/tests/test_really_cool_feature.py::test_series[int8] PASSED          [ 72%]
    mixsea/tests/test_really_cool_feature.py::test_series[int16] PASSED         [ 81%]
    mixsea/tests/test_really_cool_feature.py::test_series[int32] PASSED         [ 90%]
    mixsea/tests/test_really_cool_feature.py::test_series[int64] PASSED         [100%]

    =============================== 9 passed, 1 skipped, 1 xfailed in 0.15s ===========

Tests that we have ``parametrized`` are now accessible via the test name, for
example we could run these with ``-k int8`` to sub-select *only* those tests
which match ``int8``.


.. code-block:: shell

    (mixsea) $ pytest mixsea/tests/test_really_cool_feature.py -v -k int8
    ========================================= test session starts =====================
    platform darwin -- Python 3.8.2, pytest-5.4.1, py-1.8.1, pluggy-0.13.1 --
    plugins: datadir-1.3.1
    collected 11 items / 9 deselected / 2 selected

    mixsea/tests/test_really_cool_feature.py::test_dtypes[int8] PASSED          [ 50%]
    mixsea/tests/test_really_cool_feature.py::test_series[int8] PASSED          [100%]

    =================================== 2 passed, 9 deselected in 0.02s ===============


Running the test suite
----------------------

The tests can then be run directly inside your Git clone by typing::

    pytest mixsea

Often it is worth running only a subset of tests first around your changes before running the
entire suite. The easiest way to do this is with::

    pytest mixsea/path/to/test.py -k regex_matching_test_name

Or with one of the following constructs::

    pytest mixsea/tests/[test-module].py
    pytest mixsea/tests/[test-module].py::[TestClass]
    pytest mixsea/tests/[test-module].py::[TestClass]::[test_method]

.. _documenting.your.code:

Documenting your code
---------------------

Changes should be reflected in the release notes located in ``HISTORY.rst``.
This file contains an ongoing change log for each release.  Add an entry to this file to
document your fix, enhancement or (unavoidable) breaking change.  Make sure to include the
GitHub issue number when adding your entry (using ``:issue:`1234```, where ``1234`` is the
issue/pull request number).

If your code is an enhancement, it is most likely necessary to add usage
examples to the existing documentation.  This can be done following the section
regarding documentation :ref:`above <contributing.documentation>`.


Contributing your changes to *mixsea*
======================================

Committing your code
--------------------

Keep style fixes to a separate commit to make your pull request more readable.

Once you've made changes, you can see them by typing::

    git status

If you have created a new file, it is not being tracked by git. Add it by typing::

    git add path/to/file-to-be-added.py

Doing 'git status' again should give something like::

    # On branch shiny-new-feature
    #
    #       modified:   /relative/path/to/file-you-added.py
    #

The following defines how a commit message should be structured:

    * A subject line with `< 72` chars.
    * One blank line.
    * Optionally, a commit message body.

Please reference the relevant GitHub issues in your commit message using ``GH1234`` or
``#1234``.  Either style is fine, but the former is generally preferred.

Now you can commit your changes in your local repository::

    git commit -m

This will prompt you to type in your commit message.

Pushing your changes
--------------------

When you want your changes to appear publicly on your GitHub page, push your
forked feature branch's commits::

    git push origin shiny-new-feature

Here ``origin`` is the default name given to your remote repository on GitHub.
You can see the remote repositories::

    git remote -v

If you added the upstream repository as described above you will see something
like::

    origin      git@github.com:yourname/mixsea.git (fetch)
    origin      git@github.com:yourname/mixsea.git (push)
    upstream    git://github.com/modscripps/mixsea.git (fetch)
    upstream    git://github.com/modscripps/mixsea.git (push)

Now your code is on GitHub, but it is not yet a part of the *mixsea* project.  For that to
happen, a pull request needs to be submitted on GitHub.

Review your code
----------------

When you're ready to ask for a code review, file a pull request. Before you do, once
again make sure that you have followed all the guidelines outlined in this document
regarding code style, tests, and documentation. You should also
double check your branch changes against the branch it was based on:

#. Navigate to your repository on GitHub -- https://github.com/your-user-name/mixsea
#. Click on ``Branches``
#. Click on the ``Compare`` button for your feature branch
#. Select the ``base`` and ``compare`` branches, if necessary. This will be ``main`` and
   ``shiny-new-feature``, respectively.

Finally, make the pull request
------------------------------

If everything looks good, you are ready to make a pull request.  A pull request is how
code from a local repository becomes available to the GitHub community and can be looked
at and eventually merged into the main version.  This pull request and its associated
changes will eventually be committed to the main branch and available in the next
release.  To submit a pull request:

#. Navigate to your repository on GitHub
#. Click on the ``Pull Request`` button
#. You can then click on ``Commits`` and ``Files Changed`` to make sure everything looks
   okay one last time
#. Write a description of your changes in the ``Preview Discussion`` tab
#. Click ``Send Pull Request``.

This request then goes to the repository maintainers, and they will review
the code. If you need to make more changes, you can make them in
your branch, add them to a new commit, push them to GitHub, and the pull request
will be automatically updated.  Pushing them to GitHub again is done by::

    git push origin shiny-new-feature

This will automatically update your pull request with the latest code and restart the
:ref:`Continuous Integration <contributing.ci>` tests.


Delete your merged branch (optional)
------------------------------------

Once your feature branch is accepted into upstream, you'll probably want to get rid of
the branch. First, merge upstream main into your branch so git knows it is safe to
delete your branch::

    git fetch upstream
    git checkout main
    git merge upstream/main

Then you can do::

    git branch -d shiny-new-feature

Make sure you use a lower-case ``-d``, or else git won't warn you if your feature
branch has not actually been merged.

The branch will still exist on GitHub, so to delete it there do::

    git push origin --delete shiny-new-feature


PR checklist
------------

- **Properly comment and document your code.** See :ref:`"Documenting your code" <documenting.your.code>`.
- **Test that the documentation builds correctly** by typing ``make docs`` in the root directory or ``make html`` in the ``docs`` directory. This is not strictly necessary, but this may be easier than waiting for CI to catch a mistake. See :ref:`"Contributing to the documentation" <contributing.documentation>`.
- **Test your code**.

    - Write new tests if needed. See :ref:`"Test-driven development/code writing" <contributing.tdd>`.
    - Test the code using `Pytest <http://doc.pytest.org/en/latest/>`_. Running all tests (type ``pytest mixsea`` or ``make test`` in the root directory) takes a while, so feel free to only run the tests you think are needed based on your PR (example: ``pytest mixsea/tests/test_really_cool_feature.py``). CI will catch any failing tests.

- **Properly format your code** and verify that it passes the formatting guidelines set by `Black <https://black.readthedocs.io/en/stable/>`_ and `Flake8 <http://flake8.pycqa.org/en/latest/>`_. See :ref:`"Code formatting" <code.formatting>`. You can use `pre-commit <https://pre-commit.com/>`_ to run these automatically on each commit.

    - Run ``black mixsea`` in the root directory. This may modify some files. Confirm and commit any formatting changes.
    - Run ``flake8`` in the root directory. If this fails, it will log an error message.
    - Alternatively, the above formatting steps are combined in ``make style``. You can also type ``make style-check`` for a dry run.

- **Push your code and** `create a PR on GitHub <https://help.github.com/en/articles/creating-a-pull-request>`_.
- **Use a helpful title for your pull request** by summarizing the main contributions rather than using the latest commit message. If this addresses an `issue <https://github.com/modscripps/mixsea/issues>`_, please `reference it <https://help.github.com/en/articles/autolinked-references-and-urls>`_.

