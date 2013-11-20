Spectrum 1D
===========

Most commonly, astronomical spectra are recorded with CCDs and the resulting data product is thus a 2 dimensional
array. Through various data reduction steps (e.g. aperture extraction, IFU data reduction) these can result in one
to three dimensional data arrays. The `astropy.nddata.NDData`-class is an ideal storage format for this data, as
it provides the possibility to

The :class:`~specutils.Spectrum1D` class is one of the core classes of the ``specutils`` package. It is a child of the
:class:`~astropy.nddata.NDData` class. :class:`Spectrum1D` stores the flux and the associated units within the framework
of :class:`~astropy.nddata.NDData` and the wavelength solution is stored within a subclass of :class:`BaseSpectrum1DWCS`.

.. warning::
    `~specutils.Spectrum1D` is currently a work-in-progress, and thus it is
    possible there will be significant API changes in later versions of
    `specutils`.

There are several ways to instantiate a :class:`Spectrum1D` object::

    >>> from specutils import Spectrum1D
    >>> import numpy as np
    >>> wave = np.arange(6000, 9000)
    >>> flux = np.random.random(3000)
    >>> spec1d = Spectrum1D.from_array(wave, flux, unit='W m-2 angstrom-1 sr-1', dispersion_unit='angstrom')

One can instantiate also from astropy tables as well as files::

    >>> spec1d = Spectrum1D.from_table(spec_table, flux_column='objecta_flux')
    >>> spec1d = Spectrum1D.from_ascii('objecta.dat')

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

**Playing with WCS**:

    >>> flux = np.random.random(1000)
    >>> wave = np.sort(np.random.random(1000))
    >>> spec1d = Spectrum1D.from_array(wave, flux)
    >>> spec1d = Spectrum1D.interpolate(Spectrum1DLinearWCS(6000, dispersion_delta=1, pixel_index=np.arange(1000)))
    >>> #Now it can be written to a fits file as it is linear encoded which is representable in FITS headers.



.. automodapi:: specutils
    :no-inheritance-diagram: