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

    >>> spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49))
    >>> spec2 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.sample(49))
    >>> spec3 = spec1 + spec2
    >>> spec3 #doctest:+SKIP
    Spectrum1D([1.42980955, 0.76450583, 0.53973912, 1.12714653, 1.46747729,
            0.98485104, 1.13618017, 1.02447445, 0.68610084, 0.85083215,
            0.57521794, 1.16854341, 1.05139223, 1.44550638, 1.67533841,
            1.41277807, 0.46008942, 1.1399328 , 0.41708163, 0.61282024,
            0.40034388, 0.95204057, 0.98280167, 1.30647318, 1.43317265,
            0.71426198, 0.58622459, 1.17063336, 1.37244261, 1.06886942,
            1.54349149, 1.15019089, 0.51719866, 1.23114699, 1.16464384,
            0.90833751, 1.04018595, 1.34931354, 1.01936352, 0.39543304,
            1.22407522, 0.34658842, 1.18760707, 1.38161461, 1.05829078,
            1.57852604, 1.13365571, 0.59304282, 1.3913748 ])
    >>> spec3.wavelength
    <Quantity [ 10.,  20.,  30.,  40.,  50.,  60.,  70.,  80.,  90., 100.,
           110., 120., 130., 140., 150., 160., 170., 180., 190., 200.,
           210., 220., 230., 240., 250., 260., 270., 280., 290., 300.,
           310., 320., 330., 340., 350., 360., 370., 380., 390., 400.,
           410., 420., 430., 440., 450., 460., 470., 480., 490.] Angstrom>


Propagation of Uncertainties
----------------------------

Arithmetic operations also support the propagation of unceratinty information.

.. code-block:: python

    >>> from astropy.nddata import StdDevUncertainty

    >>> spec1 = Spectrum1D(spectral_axis=np.arange(10) * u.nm, flux=np.random.sample(10), uncertainty=StdDevUncertainty(np.random.sample(10) * 0.1))
    >>> spec2 = Spectrum1D(spectral_axis=np.arange(10) * u.nm, flux=np.random.sample(10), uncertainty=StdDevUncertainty(np.random.sample(10) * 0.1))
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

