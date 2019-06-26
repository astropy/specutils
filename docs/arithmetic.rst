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

    >>> spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49)*u.Jy)
    >>> spec2 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49)*u.Jy)
    >>> spec3 = spec1 + spec2
    >>> spec3 #doctest:+SKIP
    <Spectrum1D(flux=<Quantity [1.85874504, 0.61621214, 0.83535041, 0.79151086, 0.45719958,
               1.31989271, 0.27674835, 0.70582565, 1.32792166, 1.2417933 ,
               0.80759589, 0.30630131, 0.67328771, 1.00480227, 0.34217072,
               0.35003   , 1.2829507 , 1.47479322, 0.69394109, 1.85822987,
               0.98505397, 1.77968666, 0.86440386, 0.8547567 , 1.49315833,
               0.66677101, 1.26206613, 1.27884653, 1.42401291, 1.256702  ,
               0.97971943, 0.82923743, 1.34281845, 0.70388224, 1.19961318,
               0.98068789, 1.02384543, 0.31453759, 1.8461661 , 0.79451612,
               1.24044751, 0.48436633, 0.31850599, 0.82791702, 1.23084396,
               0.40616336, 0.60637889, 1.88098282, 0.76899619] Jy>, spectral_axis=<Quantity [ 1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11., 12., 13.,
               14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26.,
               27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39.,
               40., 41., 42., 43., 44., 45., 46., 47., 48., 49.] nm>)>


Propagation of Uncertainties
----------------------------

Arithmetic operations also support the propagation of unceratinty information.

.. code-block:: python

    >>> from astropy.nddata import StdDevUncertainty

    >>> spec1 = Spectrum1D(spectral_axis=np.arange(10) * u.nm, flux=np.random.sample(10)*u.Jy, uncertainty=StdDevUncertainty(np.random.sample(10) * 0.1))
    >>> spec2 = Spectrum1D(spectral_axis=np.arange(10) * u.nm, flux=np.random.sample(10)*u.Jy, uncertainty=StdDevUncertainty(np.random.sample(10) * 0.1))
    >>> spec1.uncertainty #doctest:+SKIP
    StdDevUncertainty([0.04386832, 0.09909487, 0.07589192, 0.0311604 ,
                   0.07973579, 0.04687858, 0.01161918, 0.06013496,
                   0.00476118, 0.06720447])
    >>> spec2.uncertainty #doctest:+SKIP
    StdDevUncertainty([0.00889175, 0.00890437, 0.05194229, 0.08794455,
                   0.09918037, 0.04815417, 0.06464564, 0.0164324 ,
                   0.04358771, 0.08260218])
    >>> spec3 = spec1 + spec2
    >>> spec3.uncertainty #doctest:+SKIP
    StdDevUncertainty([0.04476039, 0.09949412, 0.09196513, 0.09330174,
                   0.12725778, 0.06720435, 0.06568154, 0.06233969,
                   0.04384698, 0.10648737])


Reference/API
-------------
.. automodapi:: specutils.spectra.spectrum_mixin
    :no-heading:
