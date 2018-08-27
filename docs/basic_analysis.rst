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

SNR
---

Currently, specutils supports basic signal-to-noise ratio calculations.

.. code-block:: python

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from specutils.spectra import Spectrum1D
    >>> from astropy.nddata import StdDevUncertainty
    >>> from specutils.analysis import snr

    >>> uncertainty = StdDevUncertainty(0.1*np.abs(np.random.random(50))*u.Jy)
    >>> spec = Spectrum1D(spectral_axis=np.arange(50), flux=(3+np.random.randn(50))*u.Jy, uncertainty=uncertainty)
    >>> snr(spec) #doctest:+SKIP
    <Quantity 149.97247134>

Centroid
--------

Currently, specutils supports basic centroid calculations.

.. code-block:: python

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from specutils.spectra import Spectrum1D
    >>> from astropy.nddata import StdDevUncertainty
    >>> from specutils.analysis import centroid

    >>> spec = Spectrum1D(spectral_axis=np.arange(50), flux=(3+np.random.randn(50))*u.Jy)
    >>> centroid(spec) #doctest:+SKIP
    <Quantity 24.39045495 Angstrom>

And if the spectrum contains a continuum, then it should be subtracted first:
.. code-block:: python

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from specutils.spectra import Spectrum1D
    >>> from astropy.nddata import StdDevUncertainty
    >>> from specutils.fitting import continuum
    >>> from specutils.analysis import centroid

    >>> spec = Spectrum1D(spectral_axis=np.arange(50), flux=(3+np.random.randn(50))*u.Jy)
    >>> continuum_baseline = continuum(spec)
    >>> c = centroid(spec-continuum_baseline) #doctest:+SKIP

Reference/API
-------------
.. automodapi:: specutils.analysis
    :no-heading:
