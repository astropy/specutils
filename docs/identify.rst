==============================
Identifying Spectrum1D Formats
==============================

``specutils`` provides a convenience function,
`~specutils.io.registers.identify_spectrum_format`, which attempts to guess the
`~specutils.Spectrum1D` file format from the list of registered formats, and
essentially acts as a wrapper on `~astropy.io.registry.identify_format`.

This function is useful for identifying a spectrum file format without reading the
whole file with the  `~specutils.Spectrum1D.read` method.  It uses the
same identification method as ``read`` however, so it provides a convenience
of access outside of calling ``read`` without any change in underlying functionality.
It returns the best guess as to a valid format from the list of ``Formats``
as given by `~astropy.io.registry.get_formats`.

For eample, to identify a SDSS MaNGA data cube file:

.. code-block:: python

    >>> from astropy.utils.data import download_file
    >>> from specutils.io.registers import identify_spectrum_format
    >>>
    >>> url = 'https://dr17.sdss.org/sas/dr17/manga/spectro/redux/v3_1_1/8485/stack/manga-8485-1901-LOGCUBE.fits.gz'
    >>> dd = download_file(url)  # doctest: +REMOTE_DATA
    >>> identify_spectrum_format(dd)  # doctest: +REMOTE_DATA
    'MaNGA cube'

or a JWST extracted 1d spectral file:

.. code-block:: python

    >>> from specutils.io.registers import identify_spectrum_format
    >>> path = '/data/jwst/jw00626-o030_s00000_nirspec_f170lp-g235m_x1d.fits'
    >>> identify_spectrum_format(path)  # doctest: +SKIP
    'JWST x1d'
