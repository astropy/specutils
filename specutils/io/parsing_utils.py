import numpy as np
import os
import re
import urllib
import io
import contextlib

from astropy import log
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.utils.exceptions import AstropyUserWarning
import astropy.units as u
import warnings

from specutils.spectra import Spectrum


@contextlib.contextmanager
def read_fileobj_or_hdulist(*args, **kwargs):
    """ Context manager for reading a filename or file object

    Returns
    -------
    hdulist : :class:`~astropy.io.fits.HDUList`
        Provides a generator-iterator representing the open file object handle.
    """
    # Access the fileobj or filename arg
    # Do this so identify functions are useable outside of Spectrum.read context
    try:
        fileobj = args[2]
    except IndexError:
        fileobj = args[0]

    if isinstance(fileobj, fits.hdu.hdulist.HDUList):
        if fits.util.fileobj_closed(fileobj):
            hdulist = fits.open(fileobj.name, **kwargs)
        else:
            hdulist = fileobj
    elif isinstance(fileobj, io.BufferedReader):
        hdulist = fits.open(fileobj)
    else:
        hdulist = fits.open(fileobj, **kwargs)

    try:
        yield hdulist

    # Cleanup even after identifier function has thrown an exception: rewind generic file handles.
    finally:
        if not isinstance(fileobj, fits.hdu.hdulist.HDUList):
            try:
                fileobj.seek(0)
            except (AttributeError, io.UnsupportedOperation):
                hdulist.close()


def spectrum_from_column_mapping(table, column_mapping, wcs=None, verbose=False,
                                 spectrum_kwargs=None):
    """
    Given a table and a mapping of the table column names to attributes
    on the Spectrum object, parse the information into a Spectrum.

    Parameters
    ----------
    table : :class:`~astropy.table.Table`
        The table object (e.g. returned from ``Table.read('data_file')``).

    column_mapping : dict
        A dictionary describing the relation between the table columns
        and the arguments of the `~specutils.Spectrum` class, along with unit
        information. The dictionary keys should be the table column names
        while the values should be a two-tuple where the first element is the
        associated `~specutils.Spectrum` keyword argument, and the second element is the
        unit for the file column (or ``None`` to take unit from the table header)::

            column_mapping = {'FLUX': ('flux', 'Jy'),
                              'WAVE': ('spectral_axis', 'um')}

    wcs : :class:`~astropy.wcs.WCS` or :class:`~gwcs.wcs.WCS`
        WCS object passed to the Spectrum initializer.

    verbose : bool
        Print extra info.

    spectrum_kwargs : dict, optional
        Additional keyword arguments passed to the `~specutils.Spectrum` initializer,
        e.g. from `spectral_axis_metadata_from_header`.

    Returns
    -------
    :class:`~specutils.Spectrum`
        The spectrum with 'spectral_axis', 'flux' and optionally 'uncertainty'
        as identified by ``column_mapping``.
    """
    spec_kwargs = dict(spectrum_kwargs or {})

    # Associate columns of the file with the appropriate Spectrum arguments
    for col_name, (kwarg_name, cm_unit) in column_mapping.items():
        # If the table object couldn't parse any unit information,
        # fallback to the column mapper defined unit
        tab_unit = table[col_name].unit

        if tab_unit and cm_unit is not None:
            # If the table unit is defined, retrieve the quantity array for
            # the column
            kwarg_val = u.Quantity(table[col_name], tab_unit)

            # Attempt to convert the table unit to the user-defined unit.
            if verbose:
                print(f"Attempting auto-convert of table unit '{tab_unit}' to "
                      f"user-provided unit '{cm_unit}'.")

            if not isinstance(cm_unit, u.Unit):
                cm_unit = u.Unit(cm_unit)
            if cm_unit.physical_type in ('length', 'frequency', 'energy'):
                # Spectral axis column information
                kwarg_val = kwarg_val.to(cm_unit, equivalencies=u.spectral())
            elif 'spectral flux' in str(cm_unit.physical_type):
                # Flux/error column information
                kwarg_val = kwarg_val.to(cm_unit, equivalencies=u.spectral_density(1 * u.AA))
        elif tab_unit:
            # The user has provided no unit in the column mapping, so we
            # use the unit as defined in the table object.
            kwarg_val = u.Quantity(table[col_name], tab_unit)
        elif cm_unit is not None:
            # In this case, the user has defined a unit in the column mapping
            # but no unit has been defined in the table object.
            kwarg_val = u.Quantity(table[col_name], cm_unit)
        else:
            # Neither the column mapping nor the table contain unit information.
            # This may be desired e.g. for the mask or bit flag arrays.
            kwarg_val = table[col_name]

        # Transpose > 1D data to row-major format
        if kwarg_val.ndim > 1:
            kwarg_val = kwarg_val.T

        spec_kwargs.setdefault(kwarg_name, kwarg_val)

    # Ensure that the uncertainties are a subclass of NDUncertainty
    if spec_kwargs.get('uncertainty') is not None:
        spec_kwargs['uncertainty'] = StdDevUncertainty(
            spec_kwargs.get('uncertainty'))

    return Spectrum(**spec_kwargs, wcs=wcs, meta={'header': table.meta})


def generic_spectrum_from_table(table, wcs=None, spectrum_kwargs=None):
    """
    Load spectrum from an Astropy table into a Spectrum object.
    Uses the following logic to figure out which column is which:

    * Spectral axis (dispersion) is the first column with units
      compatible with ``u.spectral()`` or with length units such as 'pix'.
      Need not be present, if a valid ``wcs`` parameter is passed.

    * Flux is taken from the first column with units compatible with
      ``u.spectral_density()``, or with other likely culprits such as
      'adu' or 'cts/s'.

    * Uncertainty comes from the next column with the same units as flux.

    Parameters
    ----------
    table : :class:`~astropy.table.Table`
        Table containing a column of ``flux``, and optionally ``spectral_axis``
        and ``uncertainty`` as defined above.
    wcs : :class:`~astropy.wcs.WCS`
        A FITS WCS object. If this is present, the machinery will fall back
        and default to using the ``wcs`` to find the dispersion information.
    spectrum_kwargs : dict, optional
        Additional keyword arguments passed to the `~specutils.Spectrum` initializer,
        e.g. from `spectral_axis_metadata_from_header`.

    Returns
    -------
    :class:`~specutils.Spectrum`
        The spectrum that is represented by the data from the columns
        as automatically identified above.

    Raises
    ------
    Warns if uncertainty has zeros or negative numbers.
    Raises IOError if it can't figure out the columns.

    """
    # Local function to find the wavelength or frequency column
    def _find_spectral_axis_column(table, columns_to_search):
        """
        Figure out which column in a table holds the spectral axis (dispersion).
        Take the first column that has units compatible with u.spectral()
        equivalencies. If none meet that criterion, look for other likely
        length units such as 'pix'.
        """
        additional_valid_units = [u.Unit('pix')]
        found_column = None

        # First, search for a column with units compatible with Angstroms
        for c in columns_to_search:
            try:
                table[c].to("AA", equivalencies=u.spectral())
                found_column = c
                break
            except Exception:
                continue

        # If no success there, check for other possible length units
        if found_column is None:
            for c in columns_to_search:
                if table[c].unit in additional_valid_units:
                    found_column = c
                    break

        return found_column

    # Local function to find the flux column
    def _find_spectral_column(table, columns_to_search, spectral_axis):
        """
        Figure out which column in a table holds the fluxes or uncertainties.
        Take the first column that has units compatible with
        u.spectral_density() equivalencies. If none meet that criterion,
        look for other likely length units such as 'adu' or 'cts/s'.
        """
        additional_valid_units = [u.Unit('adu'), u.Unit('ct/s'), u.Unit('count')]
        found_column = None

        # First, search for a column with units compatible with Jansky
        for c in columns_to_search:
            try:
                # Check for multi-D flux columns
                if table[c].ndim == 1:
                    spec_ax = spectral_axis
                else:
                    # Assume leading dimension corresponds to spectral_axis
                    spec_shape = np.ones(table[c].ndim, dtype=int)
                    spec_shape[0] = -1
                    spec_ax = spectral_axis.reshape(spec_shape)
                table[c].to("Jy", equivalencies=u.spectral_density(spec_ax))
                found_column = c
                break
            except Exception:
                continue

        # If no success there, check for other possible flux units
        if found_column is None:
            for c in columns_to_search:
                if table[c].unit in additional_valid_units:
                    found_column = c
                    break

        return found_column

    # Make a copy of the column names so we can remove them as they are found
    colnames = table.colnames.copy()

    # Use the first column that has spectral unit as the dispersion axis
    spectral_axis_column = _find_spectral_axis_column(table, colnames)

    if spectral_axis_column is None and wcs is None:
        raise IOError("Could not identify column containing the wavelength, frequency or energy")
    elif wcs is not None:
        spectral_axis = None
    else:
        spectral_axis = table[spectral_axis_column].to(table[spectral_axis_column].unit)
        colnames.remove(spectral_axis_column)

    # Use the first column that has a spectral_density equivalence as the flux
    flux_column = _find_spectral_column(table, colnames, spectral_axis)
    if flux_column is None:
        raise IOError("Could not identify column containing the flux")
    flux = table[flux_column].to(table[flux_column].unit)
    colnames.remove(flux_column)
    # For > 1D data transpose to row-major format
    if flux.ndim > 1:
        flux = flux.T

    # Use the next column with the same units as flux as the uncertainty
    # Interpret it as a standard deviation and check if it has zeros or negative values
    err_column = None
    for c in colnames:
        if table[c].unit == table[flux_column].unit:
            err_column = c
            break
    if err_column is not None:
        if table[err_column].ndim > 1:
            err = table[err_column].T
        elif flux.ndim > 1:  # Repeat uncertainties over all flux columns
            err = np.tile(table[err_column], flux.shape[0], 1)
        else:
            err = table[err_column]
        err = StdDevUncertainty(err.to(err.unit))
        if np.min(table[err_column]) <= 0.:
            warnings.warn("Standard Deviation has values of 0 or less", AstropyUserWarning)
    else:
        err = None

    # Check for mask
    if 'mask' in table.colnames:
        mask = table['mask']
        if mask.ndim > 1:
            mask = mask.T
    else:
        mask = None

    # Create the Spectrum object and return it
    if wcs is not None or spectral_axis_column is not None and flux_column is not None:
        # For > 1D spectral axis transpose to row-major format and return SpectrumCollection
        spectrum = Spectrum(flux=flux, spectral_axis=spectral_axis,
                            uncertainty=err, meta={'header': table.meta}, wcs=wcs,
                            mask=mask, **(spectrum_kwargs or {}))

    return spectrum


def _fits_identify_by_name(origin, fileinp, *args,
                           pattern=r'(?i).*\.fit[s]?$', **kwargs):
    """
    Check whether input file is FITS and matches a given name pattern.
    Utility function to construct an `identifier` for Astropy I/O Registry.

    Parameters
    ----------
    fileinp : str or file-like object
        FITS file name or object (provided from name by Astropy I/O Registry).
    pattern : regex str or re.Pattern
        File name pattern to be matched.
        Note: loaders should define a pattern sufficiently specific for their
        spectrum file types to avoid ambiguous/multiple matches.
    """
    fileobj = None
    filepath = None
    if pattern is None:
        pattern = r''
    _spec_pattern = re.compile(pattern)

    if isinstance(fileinp, str):
        filepath = fileinp
        try:
            fileobj = open(filepath, mode='rb')
        except FileNotFoundError:
            # Check if path points to valid url
            try:
                fileinp = urllib.request.urlopen(filepath)
            except ValueError:
                return False
    elif fits.util.isfile(fileinp):
        fileobj = fileinp
        filepath = fileobj.name

    # Check for `urlopen` object - can only probe content if seekable
    if hasattr(fileinp, 'url') and hasattr(fileinp, 'seekable'):
        filepath = urllib.parse.unquote(fileinp.url)
        if fileinp.seekable():
            fileobj = fileinp

    check = (_spec_pattern.match(os.path.basename(filepath)) is not None and
             fits.connect.is_fits(origin, filepath, fileobj, *args))

    if fileobj is not None:
        fileobj.close()

    return check


# Header keyword pairs commonly used for the target position, in order of preference.
TARGET_KEYS = (('RA_TARG', 'DEC_TARG'), ('TARG_RA', 'TARG_DEC'), ('OBJRA', 'OBJDEC'),
               ('OBJCTRA', 'OBJCTDEC'), ('RA_OBJ', 'DEC_OBJ'), ('PLUG_RA', 'PLUG_DEC'),
               ('CAT-RA', 'CAT-DEC'), ('RA', 'DEC'))


def _parse_angle(value, unit):
    """Parse a header angle given as a number (degrees) or a sexagesimal string."""
    from astropy.coordinates import Angle

    if isinstance(value, str):
        value = value.strip()
        try:
            return float(value) * u.deg
        except ValueError:
            return Angle(value, unit=unit)
    return float(value) * u.deg


def _target_from_header(header, target_keys):
    from astropy.coordinates import SkyCoord, FK4, FK5, ICRS

    for ra_key, dec_key in target_keys:
        ra, dec = header.get(ra_key), header.get(dec_key)
        if ra is None or dec is None or ra == '' or dec == '':
            continue
        try:
            ra = _parse_angle(ra, u.hourangle)
            dec = _parse_angle(dec, u.deg)
        except (ValueError, TypeError, u.UnitsError) as err:
            log.debug(f"Could not parse target position from {ra_key}/{dec_key}: {err}")
            continue

        radesys = str(header.get('RADESYS', header.get('RADECSYS', 'ICRS'))).strip().upper()
        equinox = header.get('EQUINOX')
        if radesys.startswith('FK5'):
            frame = FK5(equinox=f'J{equinox}') if equinox else FK5()
        elif radesys.startswith('FK4'):
            frame = FK4(equinox=f'B{equinox}') if equinox else FK4()
        else:
            frame = ICRS()
        return SkyCoord(ra=ra, dec=dec, frame=frame)
    return None


def _time_from_header(header, key, scale):
    """Parse a DATE-like (ISO) or MJD-like (float) header value into a Time."""
    from astropy.time import Time

    value = header.get(key)
    if value is None or value == '':
        return None
    if isinstance(value, str):
        value = value.strip()
        if key.startswith('DATE') and 'T' not in value and len(value) == 10:
            time_key = 'TIME' + key[4:]
            if header.get(time_key):
                value = f"{value}T{str(header[time_key]).strip()}"
        return Time(value, scale=scale)
    return Time(float(value), format='mjd', scale=scale)


def _obstime_from_header(header, time_key=None):
    """
    The mid-point of the exposure: DATE-AVG/MJD-AVG, the mid-point of
    DATE-BEG/DATE-END or MJD-BEG/MJD-END, or the start (DATE-OBS, MJD-OBS or
    MJD) plus half of EXPTIME.
    """
    from astropy.time import Time

    scale = str(header.get('TIMESYS', 'UTC')).strip().lower()
    if scale not in Time.SCALES:
        scale = 'utc'

    try:
        if time_key is not None:
            return _time_from_header(header, time_key, scale)

        for key in ('DATE-AVG', 'MJD-AVG'):
            obstime = _time_from_header(header, key, scale)
            if obstime is not None:
                return obstime

        for beg_key, end_key in (('DATE-BEG', 'DATE-END'), ('MJD-BEG', 'MJD-END')):
            beg = _time_from_header(header, beg_key, scale)
            end = _time_from_header(header, end_key, scale)
            if beg is not None and end is not None:
                return beg + (end - beg) / 2

        for key in ('DATE-OBS', 'MJD-OBS', 'MJD'):
            start = _time_from_header(header, key, scale)
            if start is not None:
                exptime = header.get('EXPTIME', header.get('EXPOSURE'))
                if exptime is not None and exptime != '':
                    return start + float(exptime) * u.s / 2
                return start
    except (ValueError, TypeError) as err:
        log.debug(f"Could not parse the observation time from the header: {err}")
    return None


def _location_from_header(header, location=None):
    from astropy.coordinates import EarthLocation

    if isinstance(location, EarthLocation):
        return location
    if isinstance(location, str):
        return EarthLocation.of_site(location)
    if location is not None:
        raise TypeError("location must be an EarthLocation or the name of a site")

    xyz = [header.get(f'OBSGEO-{axis}') for axis in 'XYZ']
    if all(value is not None for value in xyz):
        return EarthLocation.from_geocentric(*[float(v) for v in xyz], unit=u.m)
    lbh = [header.get(f'OBSGEO-{axis}') for axis in 'LBH']
    if all(value is not None for value in lbh):
        return EarthLocation.from_geodetic(lon=float(lbh[0]) * u.deg, lat=float(lbh[1]) * u.deg,
                                           height=float(lbh[2]) * u.m)
    return None


def spectral_axis_metadata_from_header(header, medium=None, frame=None, location=None,
                                       target_keys=None, time_key=None):
    """
    Read the medium, reference frame, target position, observation time and
    observatory location of a spectrum from a FITS header, for use as keyword
    arguments to `~specutils.Spectrum`.

    This is intended for loaders: pass what is known about the data format
    explicitly (for instance ``medium='vacuum', frame='BARYCENT'``) and let
    the standard header keywords supply the rest.

    Parameters
    ----------
    header : `~astropy.io.fits.Header` or dict-like
        The header to read from.
    medium : `~specutils.spectra.spectral_frame.SpectralMedium`, {'vacuum', 'air'} or dict, optional
        The medium of the wavelengths. If not given, it is inferred from a
        spectral ``CTYPEn`` (or ``TCTYPn`` for a table column) when present:
        ``'AWAV'`` is air; ``'WAVE'``, ``'FREQ'``, ``'ENER'`` and ``'WAVN'``
        are vacuum.
    frame : str, optional
        FITS ``SPECSYS`` code of the reference frame. If not given, the
        ``SPECSYS`` keyword is used when present.
    location : `~astropy.coordinates.EarthLocation` or str, optional
        Where the spectrum was recorded, or the name of a site known to
        `~astropy.coordinates.EarthLocation.of_site`. If not given, the
        ``OBSGEO-[XYZ]`` or ``OBSGEO-[LBH]`` keywords are used when present.
    target_keys : tuple of str or list of tuple, optional
        The ``(RA, DEC)`` keyword pair(s) holding the target position, given as
        numbers in degrees or sexagesimal strings, in the frame given by
        ``RADESYS`` (ICRS by default). Defaults to a list of common pairs,
        see ``TARGET_KEYS``.
    time_key : str, optional
        Keyword holding the mid-point of the observation, as an ISO date/time
        string or an MJD. If not given, the mid-point is taken from
        ``DATE-AVG`` or ``MJD-AVG``, or the middle of ``DATE-BEG``/``DATE-END``
        or ``MJD-BEG``/``MJD-END``, or ``DATE-OBS`` (plus ``TIME-OBS`` if the
        date has no time), ``MJD-OBS`` or ``MJD`` plus half of ``EXPTIME``
        (or the start of the exposure if there is no ``EXPTIME``). Times are
        in the ``TIMESYS`` scale, UTC by default.

    Returns
    -------
    dict
        The ``medium``, ``frame``, ``target``, ``obstime`` and ``location``
        that could be determined; entries that could not are omitted.
    """
    from ..spectra.spectral_frame import SpectralMedium, normalize_frame
    from ..utils.wcs_utils import _CTYPE_MEDIUM

    metadata = {}

    medium = SpectralMedium.from_input(medium)
    if medium is None:
        # Image axes (CTYPEn) or table columns (TCTYPn)
        for key in [f'{prefix}{i}' for i in range(1, 10) for prefix in ('CTYPE', 'TCTYP')]:
            ctype = str(header.get(key, ''))[:4]
            if ctype in _CTYPE_MEDIUM:
                medium = SpectralMedium(_CTYPE_MEDIUM[ctype])
                break
    if medium is not None:
        metadata['medium'] = medium

    if frame is None and header.get('SPECSYS'):
        try:
            frame = normalize_frame(str(header['SPECSYS']))
        except ValueError:
            warnings.warn(f"Ignoring unrecognised SPECSYS '{header['SPECSYS']}' in header.",
                          AstropyUserWarning)
    if frame is not None:
        metadata['frame'] = normalize_frame(frame)

    if target_keys is None:
        target_keys = TARGET_KEYS
    elif isinstance(target_keys[0], str):
        target_keys = [target_keys]
    target = _target_from_header(header, target_keys)
    if target is not None:
        metadata['target'] = target

    obstime = _obstime_from_header(header, time_key=time_key)
    if obstime is not None:
        metadata['obstime'] = obstime

    location = _location_from_header(header, location=location)
    if location is not None:
        metadata['location'] = location

    return metadata
