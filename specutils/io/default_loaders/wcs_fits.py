import warnings
import _io
import logging

from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS, _wcs
from astropy.modeling import models
from astropy.utils.exceptions import AstropyUserWarning

import numpy as np
import shlex

from ...spectra import Spectrum1D, SpectrumCollection
from ..registers import data_loader, custom_writer
from ..parsing_utils import read_fileobj_or_hdulist

__all__ = ['wcs1d_fits_loader', 'non_linear_wcs1d_fits', 'non_linear_multispec_fits']

log = logging.getLogger(__name__)


def identify_wcs1d_fits(origin, *args, **kwargs):
    """
    Check whether given input is FITS and has WCS definition of WCSDIM=1 in
    specified (default primary) HDU. This is used for Astropy I/O Registry.
    On writing check if filename conforms to naming convention for this format.
    """
    whdu = kwargs.get('hdu', 1)
    # Default FITS format is BINTABLE in 1st extension HDU, unless IMAGE is
    # indicated via naming pattern or (explicitly) selecting primary HDU.
    if origin == 'write':
        return ((args[0].endswith(('wcs.fits', 'wcs1d.fits', 'wcs.fit')) or
                 (args[0].endswith(('.fits', '.fit')) and whdu == 0)) and not
                hasattr(args[2], 'uncertainty'))

    hdu = kwargs.get('hdu', 0)
    # Check if number of axes is one and dimension of WCS is one
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return (hdulist[hdu].header.get('WCSDIM', 1) == 1 and
                (hdulist[hdu].header['NAXIS'] == 1 or
                 hdulist[hdu].header.get('WCSAXES', 0) == 1 )and not
                hdulist[hdu].header.get('MSTITLE', 'n').startswith('2dF-SDSS LRG') and not
                # Check in CTYPE1 key for linear solution (single spectral axis)
                hdulist[hdu].header.get('CTYPE1', 'w').upper().startswith('MULTISPE'))


@data_loader("wcs1d-fits", identifier=identify_wcs1d_fits,
             dtype=Spectrum1D, extensions=['fits', 'fit'], priority=5)
def wcs1d_fits_loader(file_obj, spectral_axis_unit=None, flux_unit=None,
                      hdu=0, **kwargs):
    """
    Loader for single spectrum-per-HDU spectra in FITS files, with the spectral
    axis stored in the header as FITS-WCS.  The flux unit of the spectrum is
    determined by the 'BUNIT' keyword of the HDU (if present), while the
    spectral axis unit is set by the WCS's 'CUNIT'.

    Parameters
    ----------
    file_obj : str, file-like or HDUList
        FITS file name, object (provided from name by Astropy I/O Registry),
        or HDUList (as resulting from astropy.io.fits.open()).
    spectral_axis_unit : :class:`~astropy.units.Unit` or str, optional
        Units of the spectral axis. If not given (or None), the unit will be
        inferred from the CUNIT in the WCS. Note that if this is provided it
        will *override* any units the CUNIT specifies.
        The WCS CUNIT will be obtained by default from the header CUNIT1 card;
        if missing, the loader will try to extract it from the WAT1_001 card.
    flux_unit : :class:`~astropy.units.Unit` or str, optional
        Units of the flux for this spectrum. If not given (or None), the unit
        will be inferred from the BUNIT keyword in the header. Note that this
        unit will attempt to convert from BUNIT if BUNIT is present.
    hdu : int
        The index of the HDU to load into this spectrum.

    Returns
    -------
    :class:`~specutils.Spectrum1D`

    Notes
    -----
    Loader contributed by Kelle Cruz.
    """
    log.info("Spectrum file looks like wcs1d-fits")

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        header = hdulist[hdu].header
        wcs = WCS(header)

        if 'BUNIT' in header:
            data = u.Quantity(hdulist[hdu].data, unit=header['BUNIT'])
            if flux_unit is not None:
                data = data.to(flux_unit)
        else:
            data = u.Quantity(hdulist[hdu].data, unit=flux_unit)

    if spectral_axis_unit is not None:
        wcs.wcs.cunit[0] = str(spectral_axis_unit)
    elif wcs.wcs.cunit[0] == '' and 'WAT1_001' in header:
        # Try to extract from IRAF-style card or use Angstrom as default.
        wat_dict = dict((rec.split('=') for rec in header['WAT1_001'].split()))
        unit = wat_dict.get('units', 'Angstrom')
        if hasattr(u, unit):
            wcs.wcs.cunit[0] = unit
        else:  # try with unit name stripped of excess plural 's'...
            wcs.wcs.cunit[0] = unit.rstrip('s')
        log.info(f"Extracted spectral axis unit '{unit}' from 'WAT1_001'")
    elif wcs.wcs.cunit[0] == '':
        wcs.wcs.cunit[0] = 'Angstrom'

    # Compatibility attribute for lookup_table (gwcs) WCS
    wcs.unit = tuple(wcs.wcs.cunit)

    meta = {'header': header}

    if wcs.naxis > 4:
        raise ValueError('FITS file input to wcs1d_fits_loader is > 4D')
    elif wcs.naxis > 1:
        for i in range(wcs.naxis - 1, 0, -1):
            try:
                wcs = wcs.dropaxis(i)
            except(_wcs.NonseparableSubimageCoordinateSystemError) as e:
                raise ValueError(f'WCS cannot be reduced to 1D: {e} {wcs}')

    return Spectrum1D(flux=data, wcs=wcs, meta=meta)


@custom_writer("wcs1d-fits")
def wcs1d_fits_writer(spectrum, file_name, hdu=0, update_header=False, **kwargs):
    """
    Write spectrum with spectral axis defined by its WCS to (primary)
    IMAGE_HDU of a FITS file.

    Parameters
    ----------
    spectrum : :class:`~specutils.Spectrum1D`
    file_name : str
        The path to the FITS file
    hdu : int
        Header Data Unit in FITS file to write to (base 0; default primary HDU)
    update_header : bool
        Update FITS header with all compatible entries in `spectrum.meta`
    unit : str or :class:`~astropy.units.Unit`
        Unit for the flux (and associated uncertainty)
    dtype : str or :class:`~numpy.dtype`
        Floating point type for storing flux array
    """
    # Create HDU list from WCS
    try:
        wcs = spectrum.wcs
        hdulist = wcs.to_fits()
        header = hdulist[0].header
    except AttributeError as err:
        raise ValueError(f'Only Spectrum1D objects with valid WCS can be written as wcs1d: {err}')

    # Verify spectral axis constructed from WCS
    wl = spectrum.spectral_axis
    dwl = (wcs.all_pix2world(np.arange(len(wl)), 0) - wl.value) / wl.value
    if np.abs(dwl).max() > 1.e-10:
        m = np.abs(dwl).argmax()
        raise ValueError('Relative difference between WCS spectral axis and'
                         f'spectral_axis at {m:}: dwl[m]')

    if update_header:
        hdr_types = (str, int, float, complex, bool,
                     np.floating, np.integer, np.complexfloating, np.bool_)
        header.update([keyword for keyword in spectrum.meta.items() if
                       (isinstance(keyword[1], hdr_types) and
                        keyword[0] not in ('NAXIS', 'NAXIS1', 'NAXIS2'))])

    # Cannot include uncertainty in IMAGE_HDU - maybe provide option to
    # separately write this to BINARY_TBL extension later.
    if spectrum.uncertainty is not None:
        warnings.warn("Saving uncertainties in wcs1d format is not yet supported!",
                      AstropyUserWarning)

    # Add flux array and unit
    ftype = kwargs.pop('dtype', spectrum.flux.dtype)
    funit = u.Unit(kwargs.pop('unit', spectrum.flux.unit))
    flux = spectrum.flux.to(funit, equivalencies=u.spectral_density(wl))
    hdulist[0].data = flux.value.astype(ftype)

    if hasattr(funit, 'long_names') and len(funit.long_names) > 0:
        comment = f'[{funit.long_names[0]}] {funit.physical_type}'
    else:
        comment = f'[{funit.to_string()}] {funit.physical_type}'
    header.insert('CRPIX1', card=('BUNIT', f'{funit}', comment))

    # If hdu > 0 selected, prepend empty HDUs
    # Todo: implement `update` mode to write to existing files
    while len(hdulist) < hdu + 1:
        hdulist.insert(0, fits.ImageHDU())

    hdulist.writeto(file_name, **kwargs)


def identify_iraf_wcs(origin, *args, **kwargs):
    """
    IRAF WCS identifier, checking whether input is FITS and has WCS definition
    of WCSDIM=2 in specified (default primary) HDU.
    The difference to wcs1d is that this format supports 2D WCS with non-linear
    wavelength solutions. This is used for Astropy I/O Registry.
    """

    hdu = kwargs.get('hdu', 0)
    # Check if dimension of WCS is greater one.
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return ('WAT1_001' in hdulist[hdu].header and
                'WAT2_001' in hdulist[hdu].header and not
                hdulist[hdu].header.get('MSTITLE', 'n').startswith('2dF-SDSS LRG') and not
                (hdulist[hdu].header.get('TELESCOP', 't') == 'SDSS 2.5-M' and
                 hdulist[hdu].header.get('FIBERID', 0) > 0) and
                (hdulist[hdu].header.get('WCSDIM', 1) > 1 or
                 hdulist[hdu].header.get('CTYPE1', '').upper().startswith('MULTISPE')))

    # Check if dimension of WCS is greater one.
    is_wcs = ('WAT1_001' in hdulist[0].header and
              'WAT2_001' in hdulist[0].header and not
              (hdulist[0].header.get('TELESCOP', 'LEVIATHAN') == 'SDSS 2.5-M' and
               hdulist[0].header.get('FIBERID', 0) > 0) and
              (hdulist[0].header.get('WCSDIM', 1) == 2 or
               hdulist[0].header.get('CTYPE1', '').upper().startswith('MULTISPE')))

    if not isinstance(args[2], (fits.hdu.hdulist.HDUList, _io.BufferedReader)):
        hdulist.close()

    return is_wcs


@data_loader('iraf', identifier=identify_iraf_wcs, dtype=Spectrum1D, extensions=['fits'])
def non_linear_wcs1d_fits(file_obj, **kwargs):
    """Load Spectrum1D with WCS spectral axis from FITS files written by IRAF

    Parameters
    ----------

    file_obj : str, file-like or HDUList
        FITS file name, object (provided from name by Astropy I/O Registry),
        or HDUList (as resulting from astropy.io.fits.open()).

    spectral_axis_unit : :class:`~astropy.units.Unit` or str, optional
        Spectral axis unit, default is None in which case will search for it
        in the header under the keyword 'WAT1_001'.
        Note that if provided  this will *override* any units from the header.

    flux_unit : :class:`~astropy.units.Unit` or str, optional
        Flux units, default is None. If not specified will attempt to read it
        using the keyword 'BUNIT' and if this keyword does not exist it will
        assume 'ADU'.
        Note that if provided  this will *override* any units from the header.

    Returns
    -------
    :class:`~specutils.Spectrum1D`
    """
    spectral_axis, flux, meta = _read_non_linear_iraf_fits(file_obj, **kwargs)

    if spectral_axis.ndim > 1:
        log.info(f'Read spectral axis of shape {spectral_axis.shape} - '
                     'consider loading into SpectrumCollection.')
        spectral_axis = spectral_axis.flatten()

    # Check for ascending or descending values, as flattening might have joined overlapping orders.
    ds = (spectral_axis[1:] - spectral_axis[:-1]) / (spectral_axis[-1] - spectral_axis[0])
    if ds.min() < 0:
        log.warning('Non-monotonous spectral axis found, consider loading '
                        'into SpectrumCollection.')

    if flux.shape[-1] != spectral_axis.shape[-1]:
        if np.prod(flux.shape) == spectral_axis.shape[-1]:
            flux = flux.flatten()
        else:
            raise ValueError('Spectral axis and flux dimensions do not match: '
                             f'{spectral_axis.shape} != {flux.shape}!')

    return Spectrum1D(flux=flux, spectral_axis=spectral_axis, meta=meta)


@data_loader('iraf', identifier=identify_iraf_wcs, dtype=SpectrumCollection, extensions=['fits'])
def non_linear_multispec_fits(file_obj, **kwargs):
    """Load SpectrumCollection with WCS spectral axes from FITS files written by IRAF

    Loader for files containing 2D spectra, i.e. flux rows of equal length but
    different spectral axes, such as the individual orders of an Echelle spectrum.

    Parameters
    ----------

    file_obj : str, file-like or HDUList
        FITS file name, object (provided from name by Astropy I/O Registry),
        or HDUList (as resulting from astropy.io.fits.open()).

    spectral_axis_unit : :class:`~astropy.units.Unit` or str, optional
        Spectral axis unit, default is None in which case will search for it
        in the header under the keyword 'WAT1_001'.
        Note that if provided this will *override* any units from the header.

    flux_unit : :class:`~astropy.units.Unit` or str, optional
        Flux units, default is None. If not specified will attempt to read it
        using the keyword 'BUNIT' and if this keyword does not exist it will
        assume 'ADU'.
        Note that if provided this will *override* any units from the header.

    Returns
    -------
    :class:`~specutils.SpectrumCollection`
    """
    spectral_axis, flux, meta = _read_non_linear_iraf_fits(file_obj, **kwargs)

    if spectral_axis.ndim == 1:
        log.info(f'Read 1D spectral axis of length {spectral_axis.shape[0]} - '
                     'consider loading into Spectrum1D.')
        spectral_axis = spectral_axis.reshape((1, -1))
    if flux.ndim == 1:
        flux = flux.reshape((1, -1))

    if spectral_axis.shape != flux.shape:
        raise ValueError('Spectral axis and flux dimensions do not match: '
                         f'{spectral_axis.shape} != {flux.shape}!')

    return SpectrumCollection(flux=flux, spectral_axis=spectral_axis, meta=meta)


def _read_non_linear_iraf_fits(file_obj, spectral_axis_unit=None, flux_unit=None, **kwargs):
    """Read spectrum data with WCS spectral axis from FITS files written by IRAF

    IRAF does not strictly follow the fits standard especially for non-linear
     wavelength solutions.

    Parameters
    ----------

    file_obj : str, file-like or HDUList
        FITS file name, object (provided from name by Astropy I/O Registry),
        or HDUList (as resulting from astropy.io.fits.open()).

    spectral_axis_unit : :class:`~astropy.Unit` or str, optional
        Spectral axis unit, default is None in which case will search for it
        in the header under the keyword 'WAT1_001', and if none found there,
        will assume 'Angstrom'.
        Note that if provided this will *override* any units from the header.

    flux_unit : :class:`~astropy.Unit` or str, optional
        Flux units, default is None. If not specified will attempt to read it
        using the keyword 'BUNIT' and if this keyword does not exist will
        assume 'ADU'.
        Note that if provided this will *override* any units from the header.

    Returns
    -------
    Tuple of data to pass to SpectrumCollection() or Spectrum1D():

    spectral_axis : :class:`~astropy.units.Quantity`
        The spectral axis or axes as constructed from WCS(hdulist[0].header).
    flux : :class:`~astropy.units.Quantity`
        The flux data from hdulist[0].data.
    meta : dict
        Dictionary of {'header': hdulist[0].header}.
    """

    log.info('Loading 1D non-linear fits solution')

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        header = hdulist[0].header
        for wcsdim in range(1, header['WCSDIM'] + 1):
            ctypen = header['CTYPE{:d}'.format(wcsdim)]
            if ctypen == 'LINEAR':
                log.info("linear Solution: Try using `format='wcs1d-fits'` instead")
                wcs = WCS(header)
                spectral_axis = _read_linear_iraf_wcs(wcs=wcs, dc_flag=header['DC-FLAG'])
            elif ctypen == 'MULTISPE':
                log.info("Multi spectral or non-linear solution")
                spectral_axis = _read_non_linear_iraf_wcs(header=header, wcsdim=wcsdim)
            else:
                raise NotImplementedError

        if flux_unit is not None:
            data = hdulist[0].data * u.Unit(flux_unit)
        elif 'BUNIT' in header:
            data = u.Quantity(hdulist[0].data, unit=header['BUNIT'])
        else:
            log.info("Flux unit was not provided, nor found in the header. Assuming ADU.")
            data = u.Quantity(hdulist[0].data, unit='adu')

    if spectral_axis_unit is None:
        # Try to extract from IRAF-style card or use Angstrom as default.
        wat_dict = dict((rec.split('=') for rec in header['WAT1_001'].split()))
        unit = wat_dict.get('units', 'Angstrom')
        if hasattr(u, unit):
            spectral_axis_unit = unit
        else:  # try with unit name stripped of excess plural 's'...
            spectral_axis_unit = unit.rstrip('s')
        log.info(f"Extracted spectral axis unit '{spectral_axis_unit}' from 'WAT1_001'")
    spectral_axis *= u.Unit(spectral_axis_unit)

    return spectral_axis, data, dict(header=header)


def _read_linear_iraf_wcs(wcs, dc_flag):
    """Linear solution reader

    This method read the appropriate keywords. Calls the method _set_math_model
    which decides what is the appropriate mathematical model to be used and
    creates and then evaluates the model for an array.

    Parameters
    ----------

    wcs : :class:`~astropy.wcs.WCS`
        Contains wcs information extracted from the header

    dc_flag : int
        Extracted from the header under the keyword DC-FLAG which defines what
        kind of solution is described. For linear solutions it is 0 or 1.

    Returns
    -------

    spectral_axis : :class:`~numpy.ndarray`
        Mathematical model of wavelength solution evaluated for each pixel
        position

    """
    wcs_dict = {'crval': wcs.wcs.crval[0],
                'crpix': wcs.wcs.crpix[0],
                'cdelt': wcs.wcs.cd[0],
                'dtype': dc_flag,
                'pnum': wcs._naxis[0]}

    math_model = _set_math_model(wcs_dict=wcs_dict)

    spectral_axis = math_model(range(wcs_dict['pnum']))

    return spectral_axis


def _read_non_linear_iraf_wcs(header, wcsdim):
    """Read non-linear wavelength solutions written by IRAF

    Extracts the appropriate information and organize it in a dictionary for
    calling the method _set_math_model which decides what is the appropriate
    mathematical model to be used according the the type of wavelength solution
    it is dealing with.

    Parameters
    ----------

    header : :class:`~astropy.io.fits.header.Header`
        Full header of file being loaded

    wcsdim : int
        Number of the wcs dimension to be read.

    Returns
    -------

    spectral_axis : :class:`~numpy.ndarray`
        Mathematical model of wavelength solution evaluated for each pixel
        position
    """

    wat_wcs_dict = {}
    wcs_parser = {'aperture': int, 'beam': int, 'dtype': int,
                  'dstart': float, 'avdelt': float,
                  'pnum': lambda x: int(float(x)), 'z': float,
                  'alow': lambda x: int(float(x)), 'ahigh': lambda x: int(float(x)),
                  'weight': float, 'zeropoint': float,
                  'ftype': int, 'order': lambda x: int(float(x)),
                  'pmin': lambda x: int(float(x)), 'pmax': lambda x: int(float(x))}

    ctypen = header['CTYPE{:d}'.format(wcsdim)]
    log.info('Attempting to read CTYPE{:d}: {:s}'.format(wcsdim, ctypen))
    if ctypen == 'MULTISPE':
        # This is extracting all header cards for f'WAT{wcsdim}_*' into a list
        wat_head = header['WAT{:d}*'.format(wcsdim)]
        if len(wat_head) == 1:
            log.debug('Get units')
            wat_array = wat_head[0].split(' ')
            for pair in wat_array:
                split_pair = pair.split('=')
                wat_wcs_dict[split_pair[0]] = split_pair[1]
        elif len(wat_head) > 1:
            wat_string = ''
            for key in wat_head:
                wat_string += f'{header[key]:68s}'  # Keep header from stripping trailing blanks!
            wat_array = shlex.split(wat_string.replace('=', ' '))
            if len(wat_array) % 2 == 0:
                for i in range(0, len(wat_array), 2):
                    # if wat_array[i] not in wcs_dict.keys():
                    wat_wcs_dict[wat_array[i]] = wat_array[i + 1]
                    # print(wat_array[i], wat_array[i + 1])

    for key in wat_wcs_dict.keys():
        log.debug("{:d} -{:s}- {:s}".format(wcsdim, key, wat_wcs_dict[key]))

    specn = [k for k in wat_wcs_dict.keys() if k.startswith('spec')]
    spectral_axis = np.empty((len(specn), header['NAXIS1']))
    for n, sp in enumerate(specn):
        spec = wat_wcs_dict[sp].split()
        wcs_dict = dict((k, wcs_parser[k](spec[i])) for i, k in enumerate(wcs_parser.keys()))
        wcs_dict['fpar'] = [float(i) for i in spec[15:]]

        log.debug(f'Retrieving model for {sp}: {wcs_dict["dtype"]} {wcs_dict["ftype"]}')
        math_model = _set_math_model(wcs_dict=wcs_dict)

        spectral_axis[n] = math_model(range(1, wcs_dict['pnum'] + 1)) / (1. + wcs_dict['z'])

    log.info(f'Constructed spectral axis of shape {spectral_axis.shape}')
    return spectral_axis


def _set_math_model(wcs_dict):
    """Defines a mathematical model of the wavelength solution

    Uses 2 keywords to decide which model is to be built and calls the
    appropriate function.

    dtype:

        -1: None, no wavelength solution available
        0: Linear wavelength solution
        1: Log-Linear wavelength solution (not implemented)
        2: Non-Linear solutions
            ftype:
                1: Chebyshev
                2: Legendre
                3: Linear Spline (not implemented)
                4: Cubic Spline (not implemented)
                5: Pixel Coordinates (not implemented)

    Not implemented models could be implemented on user-request.


    Parameters
    ----------

    wcs_dict : dict
        Contains all the necessary wcs information needed for building any of
        models supported.


    Returns
    -------

    The mathematical model which describes the transformation from pixel to
    wavelength. An instance of `~astropy.modeling.Model`.

    """
    if wcs_dict['dtype'] == -1:
        log.debug('No wavelength solution found (DTYPE={dtype:d})'.format(**wcs_dict))
        return _none()
    elif wcs_dict['dtype'] == 0:
        log.debug('Setting model for DTYPE={dtype:d}'.format(**wcs_dict))
        return _linear_solution(wcs_dict=wcs_dict)
    elif wcs_dict['dtype'] == 1:
        log.debug('Setting model for DTYPE={dtype:d}'.format(**wcs_dict))
        return _log_linear(wcs_dict=wcs_dict)
    elif wcs_dict['dtype'] == 2:
        log.debug('Setting model for DTYPE={dtype:d} FTYPE={ftype:d}'.format(**wcs_dict))
        if wcs_dict['ftype'] == 1:
            return _chebyshev(wcs_dict=wcs_dict)
        elif wcs_dict['ftype'] == 2:
            return _non_linear_legendre(wcs_dict=wcs_dict)
        elif wcs_dict['ftype'] == 3:
            return _non_linear_cspline(wcs_dict=wcs_dict)
        elif wcs_dict['ftype'] == 4:
            return _non_linear_lspline(wcs_dict=wcs_dict)
        elif wcs_dict['ftype'] == 5:
            # pixel coordinates
            raise NotImplementedError
        elif wcs_dict['ftype'] == 6:
            # sampled coordinate array
            raise NotImplementedError
        else:
            raise SyntaxError('ftype {:d} is not defined in the '
                              'standard'.format(wcs_dict['ftype']))
    else:
        raise SyntaxError('dtype {:d} is not defined in the '
                          'standard'.format(wcs_dict['dtype']))


def _none():
    """Required to handle No-wavelength solution

    No wavelength solution is considered in the FITS standard (dtype = -1)

    This will return the identity function. It does not use
    `~astropy.modeling.models.Identity` because is not simpler to instantiate.
    Instead it uses `~astropy.modeling.models.Linear1D`

    Rretuns
    -------

        A mathematical model instance of `~astropy.modeling.models.Linear1D`
        with slope 1 and intercept 0.
    """
    model = models.Linear1D(slope=1, intercept=0)
    return model


def _linear_solution(wcs_dict):
    """Constructs a Linear1D model based on the WCS information obtained
    from the header.
    """
    intercept = wcs_dict['crval'] - (wcs_dict['crpix'] - 1) * wcs_dict['cdelt']
    model = models.Linear1D(slope=wcs_dict['cdelt'], intercept=intercept)

    return model


def _log_linear(wcs_dict):
    """Returns a log linear model of the wavelength solution.

    Not implemented

    Raises
    ------
        NotImplementedError
    """
    raise NotImplementedError


def _chebyshev(wcs_dict):
    """Returns a chebyshev model of the wavelength solution.

    Constructs a Chebyshev1D mathematical model

    Parameters
    ----------

    wcs_dict : dict
        Dictionary containing all the wcs information decoded from the header and
        necessary for constructing the Chebyshev1D model.

    Returns
    -------

        `~astropy.modeling.Model`

    """
    model = models.Chebyshev1D(degree=wcs_dict['order'] - 1,
                               domain=[wcs_dict['pmin'], wcs_dict['pmax']], )

    new_params = [wcs_dict['fpar'][i] for i in range(wcs_dict['order'])]
    model.parameters = new_params

    return model


def _non_linear_legendre(wcs_dict):
    """Returns a legendre model

    Constructs a Legendre1D mathematical model

    Parameters
    ----------

    wcs_dict : dict
        Dictionary containing all the wcs information decoded from the header and
        necessary for constructing the Legendre1D model.

    Returns
    -------

        :class:`~astropy.modeling.Model`

    """
    model = models.Legendre1D(degree=wcs_dict['order'] - 1,
                              domain=[wcs_dict['pmin'], wcs_dict['pmax']], )

    new_params = [wcs_dict['fpar'][i] for i in range(wcs_dict['order'])]
    model.parameters = new_params

    return model


def _non_linear_lspline(wcs_dict):
    """Returns a linear spline model of the wavelength solution

    Not implemented

    This function should extract certain parameters from the `wcs_dict`
    parameter and construct a mathematical model that makes the conversion from
    pixel to wavelength. All the necessary information is already contained in
    the dictionary so the only work to be done is to make the instantiation of
    the appropriate subclass of `~astropy.modeling.Model`.

    Parameters
    ----------

        wcs_dict : dict
            Contains all the WCS information decoded from an IRAF fits header.

    Raises
    ------

        NotImplementedError
    """
    raise NotImplementedError('Linear spline is not implemented')


def _non_linear_cspline(wcs_dict):
    """Returns a cubic spline model of the wavelength solution.

    This function should extract certain parameters from the `wcs_dict`
    parameter and construct a mathematical model that makes the conversion from
    pixel to wavelength. All the necessary information is already contained in
    the dictionary so the only work to be done is to make the instantiation of
    the appropriate subclass of `~astropy.modeling.Model`.

    Not implemented

     Parameters
    ----------

        wcs_dict : dict
            Contains all the WCS information decoded from an IRAF fits header.

    Raises
    ------

        NotImplementedError
    """
    raise NotImplementedError('Cubic spline is not implemented')
