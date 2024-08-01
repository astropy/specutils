import numpy as np

from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.table import Table
import astropy.units as u
from astropy.wcs import WCS

from ...spectra import Spectrum1D
from ..registers import data_loader, custom_writer
from ..parsing_utils import (generic_spectrum_from_table,
                             spectrum_from_column_mapping,
                             read_fileobj_or_hdulist)

__all__ = ['tabular_fits_loader', 'tabular_fits_writer']


def identify_tabular_fits(origin, *args, **kwargs):
    # args[0] = filename
    hdu = kwargs.get('hdu', 1)
    # Check if filename conforms to naming convention and writer has not been
    # asked to write to primary (IMAGE-only) HDU
    if origin == 'write':
        return (hdu > 0 and args[0].endswith(('.fits', '.fit')) and not
                args[0].endswith(('wcs.fits', 'wcs1d.fits', 'wcs.fit')))

    # Test if fits has extension of type BinTable and check against
    # known keys of already defined specific formats
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return (len(hdulist) > 1 and
                isinstance(hdulist[hdu], fits.BinTableHDU) and not
                (hdulist[0].header.get('TELESCOP') == 'MULTI' and
                hdulist[0].header.get('HLSPACRN') == 'MUSCLES' and
                hdulist[0].header.get('PROPOSID') == 13650) and not
                (hdulist[0].header.get('TELESCOP') == 'SDSS 2.5-M' and
                hdulist[0].header.get('FIBERID') > 0) and not
                (hdulist[0].header.get('TELESCOP') == 'HST' and
                hdulist[0].header.get('INSTRUME') in ('COS', 'STIS')) and not
                hdulist[0].header.get('TELESCOP') == 'JWST')


@data_loader("tabular-fits", identifier=identify_tabular_fits,
             dtype=Spectrum1D, extensions=['fits', 'fit'], priority=6)
def tabular_fits_loader(file_obj, column_mapping=None, hdu=1, store_data_header=False, **kwargs):
    """
    Load spectrum from a FITS file tabular extension.

    Parameters
    ----------
    file_obj : str, file-like, or :class:`~astropy.io.fits.HDUList`
            FITS file name, object (provided from name by Astropy I/O Registry),
            or HDU list (as resulting from `~astropy.io.fits.open`).
    hdu : int
        The HDU of the fits file (default: 1st extension) to read from
    store_data_header : bool
        Defaults to ``False``, which stores the primary header in ``Spectrum1D.meta['header']``.
        Set to ``True`` to instead store the header from the specified data HDU.
    column_mapping : dict
        A dictionary describing the relation between the FITS file columns
        and the arguments of the `Spectrum1D` class, along with unit
        information. The dictionary keys should be the FITS file column names
        while the values should be a two-tuple where the first element is the
        associated `Spectrum1D` keyword argument, and the second element is the
        unit for the file column (or ``None`` to take unit from the table header)::

            column_mapping = {'FLUX': ('flux', 'Jy'),
                              'WAVE': ('spectral_axis', 'um')}

    **kwargs
        Additional optional keywords passed to
        :func:`~specutils.io.parsing_utils.read_fileobj_or_hdulist`, and when
        reading from a file-like object, through to :func:`~astropy.io.fits.open`.

    Returns
    -------
    data : :class:`Spectrum1D`
        The spectrum that is represented by the data in the input table.
    """
    # Parse the wcs information. The wcs will be passed to the column finding
    # routines to search for spectral axis information in the file.
    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        wcs = WCS(hdulist[hdu].header)
        tab = Table.read(hdulist[hdu])
        if store_data_header:
            tab.meta = hdulist[hdu].header
        else:
            tab.meta = hdulist[0].header

    # Minimal checks for wcs consistency with table data -
    # assume 1D spectral axis (having shape (0, NAXIS1),
    # or alternatively compare against shape of 1st column.
    if not (wcs.naxis == 1 and wcs.array_shape[-1] == len(tab) or
            wcs.array_shape == tab[tab.colnames[0]].shape):
        wcs = None

    # If no column mapping is given, attempt to parse the file using
    # unit information
    if column_mapping is None:
        return generic_spectrum_from_table(tab, wcs=wcs)

    return spectrum_from_column_mapping(tab, column_mapping, wcs=wcs)


@custom_writer("tabular-fits")
def tabular_fits_writer(spectrum, file_name, hdu=1, update_header=False, store_data_header=False, **kwargs):
    """
    Write spectrum to BINTABLE extension of a FITS file.

    Parameters
    ----------
    spectrum : :class:`Spectrum1D`
    file_name : str, file-like or `pathlib.Path`
        File to write to. If a file object, must be opened in a writeable mode.
    hdu : int
        Header Data Unit in FITS file to write to (currently only extension HDU 1)
    update_header : bool
        Write all compatible items in ``Spectrum1D.meta`` directly to FITS header;
        this will overwrite any identically named keys from ``Spectrum1D.meta['header']``.
    store_data_header : bool
        If ``True``, store ``Spectrum1D.meta['header']`` in the header of the target data HDU
        instead of the primary header (default ``False``).
    wunit : str or `~astropy.units.Unit`
        Unit for the spectral axis (wavelength or frequency-like)
    funit : str or `~astropy.units.Unit`
        Unit for the flux (and associated uncertainty)
    wtype : str or `~numpy.dtype`
        Floating point type for storing spectral axis array
    ftype : str or `~numpy.dtype`
        Floating point type for storing flux array
    hdulist : :class:`~astropy.io.fits.HDUList`
    **kwargs
        Additional optional keywords passed to :func:`~astropy.io.fits.HDUList.writeto`.
    """
    if hdu < 1:
        raise ValueError(f'FITS does not support BINTABLE extension in HDU {hdu}.')

    header = spectrum.meta.get('header', fits.header.Header()).copy()

    # Should this warn about incompatible items or overwriting existing ones?
    if update_header:
        hdr_types = (str, int, float, complex, bool,
                     np.floating, np.integer, np.complexfloating, np.bool_)
        header.update([keyword for keyword in spectrum.meta.items() if
                       isinstance(keyword[1], hdr_types)])

    # Strip header of FITS reserved keywords
    for keyword in ['NAXIS', 'NAXIS1', 'NAXIS2']:
        header.remove(keyword, ignore_missing=True)

    # Add dispersion array and unit
    wtype = kwargs.pop('wtype', spectrum.spectral_axis.dtype)
    wunit = u.Unit(kwargs.pop('wunit', spectrum.spectral_axis.unit))
    disp = spectrum.spectral_axis.to(wunit, equivalencies=u.spectral())

    # Mapping of spectral_axis types to header TTYPE1 (no "torque/work" types!)
    dispname = str(wunit.physical_type)
    if dispname == "length":
        dispname = "wavelength"
    elif "energy" in dispname:
        dispname = "energy"

    # Add flux array and unit
    ftype = kwargs.pop('ftype', spectrum.flux.dtype)
    funit = u.Unit(kwargs.pop('funit', spectrum.flux.unit))
    flux = spectrum.flux.to(funit, equivalencies=u.spectral_density(disp))

    columns = [disp.astype(wtype), flux.astype(ftype)]
    colnames = [dispname, "flux"]

    # Include uncertainty - units to be inferred from spectrum.flux
    if spectrum.uncertainty is not None:
        try:
            unc = (
                spectrum
                .uncertainty
                .represent_as(StdDevUncertainty)
                .quantity
                .to(funit, equivalencies=u.spectral_density(disp))
            )
            columns.append(unc.astype(ftype))
            colnames.append("uncertainty")
        except RuntimeWarning:
            raise ValueError("Could not convert uncertainty to StdDevUncertainty due"
                             " to divide-by-zero error.")

    # Add mask column if present
    if spectrum.mask is not None:
        # work around https://github.com/astropy/astropy/issues/11963
        # where int8 columns are written as bool, loosing information.
        # upcast to int16 instead.
        if spectrum.mask.dtype == np.int8:
            mask = spectrum.mask.astype(np.int16)
        else:
            mask = spectrum.mask

        columns.append(mask)
        colnames.append('mask')

    # For > 1D data transpose from row-major format
    for c in range(1, len(columns)):
        if columns[c].ndim > 1:
            columns[c] = columns[c].T

    tab = Table(columns, names=colnames)
    if store_data_header:
        hdu0 = fits.PrimaryHDU()
        hdu1 = fits.BinTableHDU(data=tab, header=header)
    else:
        hdu0 = fits.PrimaryHDU(header=header)
        hdu1 = fits.BinTableHDU(data=tab)

    # This will overwrite any 'EXTNAME' previously read from a valid header; should it?
    hdu1.header.update(EXTNAME='DATA')

    hdulist = fits.HDUList([hdu0, hdu1])

    # TODO: Use output_verify options to check for valid FITS
    hdulist.writeto(file_name, **kwargs)
