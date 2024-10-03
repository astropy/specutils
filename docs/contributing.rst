.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/astropy/specutils/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug"
and "help wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

Specutils could always use more documentation, whether as part of the
official specutils docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/astropy/specutils/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome.

Get Started!
------------

Ready to contribute? Here's how to set up :ref:`specutils <specutils>` for local development.

1. Fork the :ref:`specutils <specutils>` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/specutils.git

3. Install your local copy, preferably into some sort of virtual environment using your
   preferred environment manager. For example, using conda::

    $ conda create --name specutils-dev python=3.11 pip
    $ conda activate specutils-dev
    $ cd specutils/
    $ pip install -e .'[test]'

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the tests::

    $ tox -e codestyle
    $ pytest

  The tests will run on other Python versions automatically on opening a pull request,
  but if you want to attempt to run the full test suite locally before doing so you can
  run ``tox``::

    $ tox


6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 3.10 - 3.12, and for PyPy. Check
   that all required tests passed in the Github Actions CI section at the
   bottom of your pull request.

Tips
----

To run a subset of the tests, you can call ``pytest`` with a specific file
provided as input. For example. from the base directory of a cloned
``specutils`` repository, you could run::

  $ pytest specutils/tests/test_regions.py

You can also run a specific test defined within a file using the ``-k`` flag,
for example::

  $ pytest specutils/tests/test_regions.py -k test_invert
