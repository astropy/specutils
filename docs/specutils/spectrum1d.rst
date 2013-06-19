Spectrum1D
----------

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

Once a spectrum is instantiated, one can access the `flux` and `dispersion`::

    >>> spec1d.dispersion
    array([6000, 6001, 6002, ..., 8997, 8998, 8999])
    >>> spec1d.flux
    array([ 0.34272852,  0.92834782,  0.64680224, ...,  0.03348069,
            0.10291822,  0.33614334])

.. automodapi:: specutils
    :no-inheritance-diagram: