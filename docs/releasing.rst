.. highlight:: shell

====================
Release Instructions
====================

You will need to set up a gpg key (see the `astropy docs section on this <http://docs.astropy.org/en/stable/development/releasing.html#key-signing-info>`_ for more), PyPI account, and install twine before
following these steps.

1. Ensure all of the issues slated for this release on GitHub are either closed or moved to a new milestone.
2. Pull a fresh copy of the master branch from GitHub down to your local machine.
3. Update the Changelog - Move the filled out changelog headers from unreleased to a new released section with release number and date.
   Make sure you still have empty sections for the unreleased section (Can make a new commit after this step if desired).
4. Update version number at the bottom of the setup.cfg file.
5. Make a commit with this change.
6. Tag the commit you just made (replace version numbers with your new number)::

    $ git tag -s v0.5.2 -m "tagging version 0.5.2"

7. Checkout tagged version (replace version number)::

    $ git checkout v0.5.2

8. (optional but encouraged) Run test suite locally, make sure they pass.
9. Now we do the PyPI release (steps 20,21 in the `astropy release procedures <http://docs.astropy.org/en/stable/development/releasing.html>`_)::

    $ git clean -dfx
    $ cd astropy_helpers; git clean -dfx; cd ..
    $ python setup.py build sdist
    $ gpg --detach-sign -a dist/specutils-0.5.1.tar.gz
    $ twine upload dist/specutils-0.5.1.tar.gz

10. Checkout master.
11. Back to development - update setup.cfg version number back to dev, i.e. 0.6.dev and make a commit.
12. Push to Github with  “--tags” parameter (you may need to lift direct master push restrictions on the GitHub repo)
13. Do "release" with new tag on GitHub repo.
14. If there is a milestone for this release, "close" the milestone on GitHub.
15. Double-check (and fix if necessary) that relevant conda builds have proceeded sucessfully (e.g. https://github.com/conda-forge/specutils-feedstock)
