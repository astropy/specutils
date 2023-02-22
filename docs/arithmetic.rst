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
    <Spectrum1D(flux=<Quantity [0.85594057, 0.59914105, 0.86545315, 1.29308365, 0.56743587,
               0.63720232, 1.03919556, 0.33693653, 0.89068491, 1.41613598,
               0.72461457, 1.20411351, 0.96480272, 0.37496506, 0.70241888,
               1.36924151, 0.90943254, 0.82210346, 0.98018949, 1.05861761,
               0.26172516, 1.02205189, 0.51839963, 1.21572449, 0.87754143,
               1.02493144, 0.95316681, 0.37872965, 0.17723648, 1.21662474,
               1.39171024, 1.23614795, 1.10636247, 0.97294585, 1.5453743 ,
               1.01020945, 1.42125961, 1.36636734, 1.11338214, 0.58687869,
               0.63074156, 0.67475105, 0.54093389, 1.77345469, 1.22990398,
               0.78162068, 0.63760289, 0.63139356, 0.97112644] Jy>, spectral_axis=<SpectralAxis
       (observer to target:
          radial_velocity=0.0 km / s
          redshift=0.0)
      [ 1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11., 12., 13., 14.,
       15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28.,
       29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42.,
       43., 44., 45., 46., 47., 48., 49.] nm>)>


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
