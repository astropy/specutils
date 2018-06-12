==============
Basic Analysis
==============

The specutils package comes with a few basic spectral analytic functions.
More extensive and involve analysis techniques will be available in another
package, `specreduce <https://github.com/astropy/specreduce>`_.

Specutils supports some built-in callable functions for basic calculations
over the given spectrum object.

Equivalent Width
----------------

Currently, specutils supports basic equivalent width calculations.

.. code-block:: python

    >>> import numpy as np
    >>> from specutils.spectra import Spectrum1D
    >>> from specutils.analysis import equivalent_width

    >>> spec = Spectrum1D(spectral_axis=np.arange(50), flux=np.random.randn(50))
    >>> equivalent_width(spec) #doctest:+SKIP
    <Quantity 24.16006697 Angstrom>


Reference/API
-------------
.. automodapi:: specutils.analysis
    :no-heading:
