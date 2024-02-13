===================
Spectrum Arithmetic
===================

Specutils sports the ability to perform arithmetic operations over spectrum
data objects. There is full support for propagating unit information.

.. note:: Spectrum arithmetic requires that the two spectrum objects have
          compatible WCS information.

.. warning:: Specutils does not currently implement interpolation techniques
             for converting spectral axes information from one WCS source to
             another.


Basic Arithmetic
----------------

Arithmetic support includes addition, subtract, multiplication, and division.

.. code-block:: python

    >>> from specutils import Spectrum1D
    >>> import astropy.units as u
    >>> import numpy as np

    >>> rng = np.random.default_rng(12345)
    >>> spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=rng.random(49)*u.Jy)
    >>> spec2 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=rng.random(49)*u.Jy)
    >>> spec3 = spec1 + spec2
    >>> spec3  # doctest: +FLOAT_CMP
    <Spectrum1D(flux=[0.8559405665668484 ... 0.9711264429515736] Jy (shape=(49,), mean=0.91592 Jy); spectral_axis=<SpectralAxis
       (observer to target:
          radial_velocity=0.0 km / s
          redshift=0.0)
      [ 1.  2.  3. ... 47. 48. 49.] nm> (length=49))>


Propagation of Uncertainties
----------------------------

Arithmetic operations also support the propagation of unceratinty information.

.. code-block:: python

    >>> from astropy.nddata import StdDevUncertainty

    >>> rng = np.random.default_rng(12345)
    >>> wave = np.arange(10) * u.nm
    >>> spec1 = Spectrum1D(spectral_axis=wave,
    ...                    flux=rng.random(10)*u.Jy,
    ...                    uncertainty=StdDevUncertainty(rng.random(10) * 0.1))
    >>> spec2 = Spectrum1D(spectral_axis=wave,
    ...                    flux=rng.random(10)*u.Jy,
    ...                    uncertainty=StdDevUncertainty(rng.random(10) * 0.1))
    >>> spec1.uncertainty
    StdDevUncertainty([0.02482457, 0.09488812, 0.06672375, 0.00958979,
                       0.04418397, 0.08864799, 0.06974535, 0.03264729,
                       0.07339282, 0.0220135 ])
    >>> spec2.uncertainty
    StdDevUncertainty([0.08547419, 0.06016212, 0.09319884, 0.07247814,
                       0.08605513, 0.09293378, 0.0546186 , 0.0937673 ,
                       0.04949879, 0.02737732])
    >>> spec3 = spec1 + spec2
    >>> spec3.uncertainty
    StdDevUncertainty([0.08900616, 0.11235317, 0.11462147, 0.07310981,
                       0.09673525, 0.12843346, 0.08858671, 0.09928822,
                       0.08852478, 0.03512992])

Reference/API
-------------
.. automodapi:: specutils.spectra.spectrum_mixin
    :no-heading:
