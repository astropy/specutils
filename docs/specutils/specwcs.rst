Spectral World Coordinate Systems
=================================

.. currentmodule:: specutils.wcs.specwcs

In its simplest representation a spectrum is a one-dimensional array storing the flux information. In this form pixel indices
form a coordinate system for the spectrum. In most cases, however, the data has been calibrated in a way that
allows a mapping between pixel coordinates and physical coordinates. These transformations are often referred to
as a world coordinate system (WCS) and can be thought of as functions:

.. math::
    f_{\textrm{WCS}}(c_0, c_1, \dots, c_l;\ x_1, x_2, \dots, x_n) \mapsto y_1, y_2, \dots, y_m,

where :math:`c_i` represents a parameter for the transform, :math:`x_j` are pixel coordinates and :math:`y_k` are spectral coordinates (e.g. angstrom).

:class:`~astropy.modeling.core.Model` is essentially a parametrized function and thus is an ideal way to represent WCS transformations.
In addition, there are fitters available (link) that can be used to create WCS transforms from calibration data.
Here's how a simple model (:math:`f(c_0, c_1;\ x_0) = c_0 \exp{(c_1 x_0)}`) would look in python::

    class SimpleModel(astropy.modeling.Model):

        c_0 = astropy.modeling.Parameter('c_0')
        c_1 = astropy.modeling.Parameter('c_1')

        def __init__(self, c_0, c_1):
            super(SimpleModel, self).__init__()
            self.c_0 = c_0
            self.c_1 = c_1

        def __call__(self, x_0):
            y_0 = self.c_0 * np.exp(self.c_1 * x_0)
            return y_0

It would be used in the following way::

    >>> mymodel = SimpleModel(5, 6)
    >>> mymodel(7)
    8.6963747076025047e+18



Spectral 1D WCS
---------------

The simplest WCS for 1D is a lookup table. This is often found in a ascii tables where one column
is the lookup table (wavelength array) and one column is the flux array. In terms of the functional transform view:
the lookup table represents a parameter (:math:`c_0`)::

    >>> from specutils import Spectrum1D
    >>> from specutils.wcs import specwcs
    >>> import numpy as np
    >>> wave = np.arange(6000, 9000) # wavelength
    >>> flux = np.random.random(3000) # flux
    >>> lookup_table_wcs = specwcs.Spectrum1DLookupWCS(wave, unit='angstrom') # Creating the wcs
    >>> lookup_table_wcs(0) # doing the transformation for pixel 0
    <Quantity 6000.0 Angstrom>
    >>> Spectrum1D(flux=flux, wcs=lookup_table_wcs)
    Spectrum1D([ 0.66118716,  0.39584688,  0.81210479, ...,  0.5238203 ,
             0.05986459,  0.11621466])

For more information about the Spectrum1d object go to :doc:`spectrum1d`.

Another common WCS is the linear dispersion (:math:`f_{\textrm{WCS}}(x) = c_0 + c_1 x`) and commonly serialized (encoded) to FITS keyword headers.
For linear dispersion we are using the general `Spectrum1DPolynomialWCS` WCS::

    >>> from specutils.wcs import specwcs
    >>> linear_wcs = specwcs.Spectrum1DPolynomialWCS(degree=1, c0=6000, c1=1., unit='angstrom')
    >>> linear_wcs(1.5)
    <Quantity 6001.5 Angstrom>

There are also a couple of models which can be useful to combine different WCS. These are the `CompositeWCS` and
`WeightedCombinationWCS` models.

When there is a need to apply more than one WCS in a particular order to the input, `CompositeWCS` should be used::

    >>> import math
    >>> from specutils.wcs import specwcs
    >>> linear_wcs = specwcs.Spectrum1DPolynomialWCS(degree=1, c0=6000, c1=1.)
    >>> log = lambda x: math.log(x, 10)
    >>> composite_wcs = specwcs.CompositeWCS([linear_wcs, log])
    >>> composite_wcs(1.5)
    3.778259810434678
    >>> math.log(linear_wcs(1.5), 10)
    3.778259810434678
    >>> add_one = lambda x: x+1
    >>> composite_wcs.add_WCS(add_one)
    >>> composite_wcs(1.5)
    4.778259810434678

The above code applies two WCS to the input. First, it applies the polynomial WCS, and then returns the logarithm
of the output. As illustrated, one can also add user defined functions to this composite WCS. The `add_WCS` method
can also be used apart from the init function to add more WCS to the composite WCS. In general, the `CompositeWCS`
returns the following:

.. math::
    wcs_n(wcs_n-1\dots(wcs_2(wcs_1(input))\dots))

The `WeightedCombinationWCS` model provides another way of combines different WCS. In this model, each WCS gets the
same input. The output from each WCS is added together using a weight and a zero point offset::

    >>> import math
    >>> from specutils.wcs import specwcs
    >>> linear_wcs1 = specwcs.Spectrum1DPolynomialWCS(degree=1, c0=6000, c1=1.)
    >>> linear_wcs2 = specwcs.Spectrum1DPolynomialWCS(degree=1, c0=1000, c1=-1.)
    >>> combination_wcs = specwcs.WeightedCombinationWCS()
    >>> combination_wcs.add_WCS(linear_wcs1, weight=0.5)
    >>> combination_wcs.add_WCS(linear_wcs2, weight=1.2, zero_point_offset=1000)
    >>> combination_wcs([1.5, 2])
    array([ 5398.95,  5398.6 ])
    >>> 0.5 * linear_wcs1([1.5, 2]) + 1.2 * (1000 + linear_wcs2([1.5, 2]))
    array([ 5398.95,  5398.6 ])

The order of WCS does not matter in this case. The weight and the zero point offset are the parameters stored with
each WCS, and are evaluated as shown in the above example. In general the `WeightedCombinationWCS` returns the
following:

.. math::
    \sum\limits_{i=1}^n weight_i * (zero\_point\_offset_i + WCS_i(input))

Another important model available is the `DopplerShift` model. This model is specifically for calculating the doppler
shift from velocity (v). The doppler factor is computed using the following formula:

.. math::
    doppler\_factor = \sqrt{\frac{1 + \frac{v}{c}}{1 - \frac{v}{c}}}

,where c is the speed of light

When the model is called on an input, input * doppler_factor is returned. The inverse of this model can also be
obtained, which divides the input by the doppler factor. The Doppler shift model can also be instantiated from
redshift(z) and the Doppler factor itself.

For more information about reading and writing WCS transformations to the FITS File format see :doc:`FITS reader <read_fits>`
and :doc:`FITS writer <write_fits>`.

In addition, the following WCS models exist as well:

 * `~specutils.wcs.specwcs.Spectrum1DIRAFLegendreWCS`
 * `~specutils.wcs.specwcs.Spectrum1DIRAFChebyshevWCS`
 * `~specutils.wcs.specwcs.Spectrum1DIRAFBSplineWCS`
 * `~specutils.wcs.specwcs.MultispecIRAFCompositeWCS`

 These models are just IRAF extensions of already existing models. This extensions enable the user to read and write
 from IRAF FITS files. 




.. Interpolation
    -------------
    In the following, we explain interpolation in one-dimension, however the concepts remain the same in n-dimensional space.
    Interpolation means going from the current WCS transform axis(:math:`f_\textrm{WCS}`) to another WCS transform axis (:math:`g_textrm{WCS}`).
    In the first instance, one needs to know what pixel in the new WCS transform axis represents what pixel (and that might often be a fractional one) in the old one:
    For example a simple shift by half a pixel:
    .. math::
        \textrm{pixel}_\textrm{new} 0 \mapsto \textrm{pixel}_\textrm{old} 0.5; \dots
        \textrm{pixel}_\textrm{new} 1 \mapsto \textrm{pixel}_\textrm{old} 1.5\\
        \Rightarrow h(\textrm{pixel}_\textrm{new}) \mapsto \textrm{pixel}_\textrm{old}
    Knowing both transformations (:math:`f_\textrm{WCS}` and :math:`g_\textrm{WCS}`), one can create the function
    :math:`h = g(f^{-1})` and transform represent the new pixel coordinates in the old system and interpolate.
    .. warning::
        This does not currently work and is being developed
    This example showcases the API. Currently the interpolate method accepts a WCS object or an `~numpy.ndarray`
    (then constructs a `~specutils.wcs.specwcs.Spectrum1DLookupWCS` out of it)::
        >>> from specutils import Spectrum1D
        >>> from specutils.wcs import specwcs
        >>> import numpy as np
        >>> wave = np.arange(6000, 9000) # wavelength
        >>> flux = np.random.random(3000) # flux
        >>> myspec = Spectrum1D.from_array(wave, flux, dispersion_unit='angstrom')
        >>> myspec.interpolate(wave + 0.5)




Reference/API
-------------

.. automodapi:: specutils.wcs.specwcs
    :no-inheritance-diagram:
    :skip: Model, Parameter, deprecated, check_valid_unit

