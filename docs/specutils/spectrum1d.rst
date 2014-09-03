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
    >>> flux = np.random.random(3000) * u.Unit('W m-2 angstrom-1 sr-1')
    >>> spec1d = Spectrum1D.from_array(wave, flux)
    >>> spec1d.wavelength
    <Quantity [ 6000., 6001., 6002.,...,  8997., 8998., 8999.] Angstrom>
    >>> spec1d.flux
    <Quantity [ 0.75639906, 0.23677036, 0.08408417,...,  0.82740303,
            0.38345114, 0.77815595] W / (Angstrom m2 sr)>

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

A spectrum1D object can also be sliced using `slice_index` method. This method accepts the start, stop and step for
a slice. This slice works similar to the python-like slicing for lists::

    >>> alternate_spec1d = spec1d.slice_index(step=2)
    >>> alternate_spec1d.dispersion
    array([ 6000.,  6002.,  6004., ...,  8994.,  8996.,  8998.])
    >>> alternate_spec1d.flux
    array([ 0.53091768,  0.81469917,  0.85442559, ...,  0.74203082, 0.453299,  0.21887765])
    >>> reverse_spec1d = spec1d.slice_index(step=-1)
    >>> reverse_spec1d.dispersion
    array([ 8999.,  8998.,  8997., ...,  6002.,  6001.,  6000.])
    >>> mixed_spec1d = reverse_spec1d.slice_index(start=101, stop=2500, step=2)
    >>> mixed_spec1d.dispersion
    array([ 8898.,  8896.,  8894., ...,  6504.,  6502.,  6500.])
    
If not provided, this method makes the same assumptions about start, stop and step as slicing for python lists. `slice_index` is used to
slice a spectrum from its indices. The start, stop and step represent indices at which dispersion, flux and other quantities are computed.
In future, it might be possible to slice at a particular dispersion or flux value, thus the distinction.

To learn more about the WCS transformations please go to :doc:`the WCS documentation page <specwcs>`. In many cases, spectra are stored in FITS files.
FITS readers and writers are documented at :doc:`FITS reader docs<read_fits>` and :doc:`FITS writer docs<write_fits>`..

.. automodapi:: specutils
    :no-inheritance-diagram: