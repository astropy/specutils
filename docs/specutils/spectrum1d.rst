Spectrum 1D
===========

.. warning::
    `~specutils.Spectrum1D` is currently a work-in-progress, and thus it is
    possible there will be significant API changes in later versions of
    `specutils`.


The :class:`~specutils.Spectrum1D` class is one of the core classes of the ``specutils`` package. It inherits most of
its properties from the :class:`~astropy.nddata.NDData` class. :class:`Spectrum1D` stores the flux and the associated units within the framework
of :class:`~astropy.nddata.NDData` and the wavelength solution is stored within a subclass of :class:`BaseSpectrum1DWCS`.


There are several ways to instantiate a :class:`Spectrum1D` object::

    >>> from specutils import Spectrum1D
    >>> from astropy import units as u
    >>> import numpy as np
    >>> wave = np.arange(6000, 9000) * u.Angstrom
    >>> flux = np.random.random(3000) # due to current limitations in NDData, a quantity can not be given
    >>> spec1d = Spectrum1D.from_array(wave, flux, unit='W m-2 angstrom-1 sr-1')


The ability to instantiate from any two arrays allows to load from ascii files (see `astropy.io.ascii <http://docs.astropy.org/en/stable/io/ascii/index.html>`_)
and other data structures like `~numpy.recarray` and `~astropy.table.Table`.

Using linear dispersion is another very common way to store dispersion
(:math:`f(\textrm{pixel}) = c_0 + c_1 \times \textrm{pixel}`). Here a simple example::

    >>> from specutils import Spectrum1D
    >>> from specutils.wcs import specwcs
    >>> linear_wcs = specwcs.Spectrum1DPolynomialWCS(degree=1, unit='angstrom', c0=6000, c1=1.0)
    >>> Spectrum1D(flux=flux, wcs=linear_wcs)
    Spectrum1D([ 0.66118716,  0.39584688,  0.81210479, ...,  0.5238203 ,
             0.05986459,  0.11621466])

Once a spectrum is instantiated, one can access the `flux`, `wavelength`, `frequency`, `energy`
(default units can be changed with `energy_unit`, etc.)::

    >>> spec1d.wavelength
    <Quantity [ 6000., 6001., 6002.,...,  8997., 8998., 8999.] Angstrom>
    >>> spec1d.frequency
    <Quantity [  4.99654097e+14,  4.99570835e+14,  4.99487601e+14,...,
                 3.33213802e+14,  3.33176770e+14,  3.33139747e+14] Hz>
    >>> spec1d.energy
    <Quantity [  3.31074281e-19,  3.31019111e-19,  3.30963959e-19,...,
             2.20789784e-19,  2.20765246e-19,  2.20740714e-19] J>
    >>> spec1d.energy_unit = 'eV'
    >>> spec1d.energy
    <Quantity [ 2.06640322, 2.06605887, 2.06571464,...,  1.3780615 ,
            1.37790835, 1.37775523] eV>
    >>> spec1d.flux
    array([ 0.34272852,  0.92834782,  0.64680224, ...,  0.03348069,
            0.10291822,  0.33614334])

To learn more about the WCS transforms please go to :doc:`the WCS documentation page <specwcs>`. In many cases, spectra are stored in FITS files.
FITS readers are documented at :doc:`the FITS WCS documentation page <fits_wcs>`.

.. automodapi:: specutils
    :no-inheritance-diagram: