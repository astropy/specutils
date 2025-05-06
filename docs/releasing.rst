.. highlight:: shell

====================
Release Instructions
====================


1. Ensure all of the issues slated for this release on GitHub are either closed or moved to a new milestone.
2. Pull a fresh copy of the main branch from GitHub down to your local machine.
3. Update the Changelog - Move the filled out changelog headers from unreleased to a new released section with release version number.
4. Make a commit with this change and push to main, if you have permission, otherwise open a PR to main.
5. Go to `Releases on GitHub <https://github.com/spacetelescope/jdaviz/releases>`_
   and `create a new GitHub release <https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository>`_
   targeting ``main``, and give it a new ``vX.Y.Z`` tag (do not choose any existing tags).
   Copy the relevant section from CHANGES.rst into the release notes section and clean up
   any formatting problems.
6. Do "release" with new tag on GitHub repo.
7. If there is a milestone for this release, "close" the milestone on GitHub.
8. Check that the release was successfully uploaded to PyPi.
9. Double-check (and fix if necessary) that relevant conda builds have proceeded sucessfully (e.g. https://github.com/conda-forge/specutils-feedstock)
