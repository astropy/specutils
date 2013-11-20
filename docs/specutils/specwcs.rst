Spectral World Coordinate Systems
=================================

In its simplest form a spectrum is a one-dimensional array storing the flux information. In this form pixel indices
form a coordinate system for the spectrum. In most cases, however, the data has been calibrated in a way that
allows a mapping between pixel coordinates and physical coordinates. These transformations are often referred to
 as a world coordinate system (WCS) and can be thought of as functions:

.. math::
    f_{\textrm{WCS}}(x_1, x_2, \dots, x_n) \mapsto y_1, y_2, \dots, y_m,

where :math:`x_i` are pixel coordinates and :math:`y_j` are spectral coordinates (e.g. angstrom).

as function using modeling.Models. Modeling has parameters, and input and output dimensions making it perfect
for what we're looking for. In addition, fitters available. Here's an example of a simple model::





Spectral 1D WCS
---------------

The simplest (and most commonly used) WCS for 1D is a lookup table. This is often found in a ascii tables where one column
is the lookup table (wavelength array) and one column is the flux array::

    >>> from specutils import Spectrum1D
    >>> from specutils.wcs import specwcs
    >>> import numpy as np
    >>> wave = np.arange(6000, 9000) # wavelength
    >>> flux = np.random.random(3000) # flux
    >>> lookup_table_wcs = specwcs.Spectrum1DLookupWCS(wave, unit='angstrom') # Creating the wcs
    >>> lookup_table_wcs(0) # doing the transformation for pixel 0
    <Quantity 6000.0 Angstrom>
    >>> Spectrum1D(data=flux, wcs=lookup_table_wcs)
    Spectrum1D([ 0.66118716,  0.39584688,  0.81210479, ...,  0.5238203 ,
             0.05986459,  0.11621466])

For more information about the Spectrum1d object go to :doc:`spectrum1d`.

Another common WCS is the linear dispersion (:math:`f_{\textrm{WCS}(x) = c_0 + c_1 x`) and commonly serialized (encoded) to FITS keyword headers.
For linear dispersion we are using the general `Spectrum1DPolynomialWCS` WCS::

    >>> from specutils.wcs import specwcs
    >>> linear_wcs = specwcs.Spectrum1DPolynomialWCS(degree=1, c0=6000, c1=1., unit='angstrom')
    >>> linear_wcs(1.5)
    <Quantity 6001.5 Angstrom>


For more information about reading and writing WCS transformations to the FITS File format see `FITS and WCS <fitswcs>`.

In addition, the following WCS models exist as well:





Interpolation
-------------


Reference/API
-------------

.. automodapi:: specutils.wcs.specwcs
    :no-inheritance-diagram:
    :skip: Model, Parameter, deprecated, check_valid_unit

