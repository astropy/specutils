.. highlight:: shell

====================
Release Instructions
====================

You will need to set up a gpg key (see the `astropy docs section on this <http://docs.astropy.org/en/stable/development/releasing.html#key-signing-info>`_ for more), PyPI account, and install twine before
following these steps.

1. Ensure all of the issues slated for this release on GitHub are either closed or moved to a new milestone.
2. Pull a fresh copy of the main branch from GitHub down to your local machine.
3. Update the Changelog - Move the filled out changelog headers from unreleased to a new released section with release version number.
4. Make a commit with this change.
5. Tag the commit you just made (replace version numbers with your new number)::

    $ git tag -s v0.5.2 -m "tagging version 0.5.2"

6. Checkout tagged version (replace version number)::

    $ git checkout v0.5.2

7. (optional but encouraged) Run test suite locally, make sure they pass.
8. Now we do the PyPI release (steps 20,21 in the `astropy release procedures <http://docs.astropy.org/en/stable/development/releasing.html>`_)::

    $ git clean -dfx
    $ cd astropy_helpers; git clean -dfx; cd ..
    $ python setup.py build sdist
    $ gpg --detach-sign -a dist/specutils-0.5.1.tar.gz
    $ twine upload dist/specutils-0.5.1.tar.gz

9. Checkout main.
10. Back to development - add the next version number to the changelog as an "unreleased" section
11. Push to Github with  ``--tags`` (you may need to lift direct main push restrictions on the GitHub repo)
12. Do "release" with new tag on GitHub repo.
13. If there is a milestone for this release, "close" the milestone on GitHub.
14. Double-check (and fix if necessary) that relevant conda builds have proceeded sucessfully (e.g. https://github.com/conda-forge/specutils-feedstock)

