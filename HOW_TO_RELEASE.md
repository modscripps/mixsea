How to issue a mixsea release
-----------------------------
Copied and adjusted from xarray.

1. Ensure your main branch is synced to upstream:
     ```
     git pull upstream main
     ```
2. Look over HISTORY.rst and make sure it is complete for the current version.
    (check the date!) Consider adding a brief summary note describing the
    release at the top.
    Things to watch out for:
    - Important new features should be highlighted towards the top.
    - Function/method references should include links to the API docs.
    - Sometimes notes get added in the wrong section of the history, typically
      due to a bad merge. Check for these before a release by using git diff,
      e.g., `git diff v0.X.Y HISTORY.rst` where 0.X.Y is the previous
      release.
3. If you have any doubts, run the full test suite one final time!
     ```
     pytest
     ```
4. Check that the ReadTheDocs build is passing.
5. Bump the version by running
     ```
     bump2version patch # possible: major / minor / patch
     ```
     This will change the version number in `setup.py`, `setup.cfg` and `mixsea/__init.py__`.
6. On the main branch, commit the release in git:
     ```
     git commit -am 'Release v0.X.Y'
     ```
7. Tag the release:
     ```
     git tag -a v0.X.Y -m 'v0.X.Y'
     ```
8. Build source and binary wheels for pypi:
     ```
     make clean
     git clean -xdf  # this deletes all uncommited changes!
     python setup.py bdist_wheel sdist
     ```
9. Use twine to check the package build:
     ```
     twine check dist/mixsea-0.X.Y*
     ```
10. Use twine to register and upload the release on test.pypi.org:
     ```
     twine upload --repository-url https://test.pypi.org/legacy/ dist/mixsea-0.X.Y*
     ```
    You will need to be listed as a package owner at
    https://test.pypi.python.org/pypi/mixsea for this to work. 

11. Use twine to register and upload the release on pypi. Be careful, you can't
    take this back!
     ```
     twine upload dist/mixsea-0.X.Y*
     ```
    You will need to be listed as a package owner at
    https://pypi.python.org/pypi/mixsea for this to work.
12. Push your changes to main:
     ```
     git push upstream main
     git push upstream --tags
     ```
13. Update the stable branch (used by ReadTheDocs) and switch back to main:
     ```
     git checkout stable
     git rebase main
     git push upstream stable
     git checkout main
     ```
    It's OK to force push to 'stable' if necessary. (We also update the stable 
    branch with `git cherrypick` for documentation only fixes that apply the 
    current released version.)
14. Add a section for the next release (v.X.Y+1) to HISTORY.rst:
     ```
     .. _whats-new.0.X.Y+1:

     v0.X.Y+1 (unreleased)
     ---------------------

     Breaking changes
     ~~~~~~~~~~~~~~~~


     New Features
     ~~~~~~~~~~~~


     Bug fixes
     ~~~~~~~~~


     Documentation
     ~~~~~~~~~~~~~


     Internal Changes
     ~~~~~~~~~~~~~~~~
     ```
15. Commit your changes and push to main again:
      ```
      git commit -am 'New whatsnew section'
      git push upstream main
      ```
    You're done pushing to main!
16. Issue the release on GitHub. Click on "Draft a new release" at
    https://github.com/modscripps/mixsea/releases. Type in the version number, but
    don't bother to describe it -- we maintain that on the docs instead.
17. Update the docs. Login to https://readthedocs.org/projects/mixsea/versions/
    and switch your new release tag (at the bottom) from "Inactive" to "Active".
    It should now build automatically.

Note on version numbering:

We follow a rough approximation of semantic version. Only major releases (0.X.0)
should include breaking changes. Minor releases (0.X.Y) are for bug fixes and
backwards compatible new features, but if a sufficient number of new features
have arrived we will issue a major release even if there are no compatibility
breaks.

Once the project reaches a sufficient level of maturity for a 1.0.0 release, we
intend to follow semantic versioning more strictly.
