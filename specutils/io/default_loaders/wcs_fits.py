import logging
import warnings
import _io

from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS, _wcs
from astropy.modeling import models
from astropy.utils.exceptions import AstropyUserWarning

import numpy as np
import shlex

from ...spectra import Spectrum1D
from ..registers import data_loader, custom_writer
from ..parsing_utils import read_fileobj_or_hdulist

__all__ = ['wcs1d_fits_loader', 'wcs1d_fits_writer', 'non_linear_wcs1d_fits']


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
    file_obj: str, file-like or HDUList
        FITS file name, object (provided from name by Astropy I/O Registry),
        or HDUList (as resulting from astropy.io.fits.open()).
    spectral_axis_unit: str or `~astropy.Unit`, optional
        Units of the spectral axis. If not given (or None), the unit will be
        inferred from the CUNIT in the WCS. Note that if this is provided it
        will *override* any units the CUNIT provides.
    flux_unit: str or `~astropy.Unit`, optional
        Units of the flux for this spectrum. If not given (or None), the unit
        will be inferred from the BUNIT keyword in the header. Note that this
        unit will attempt to convert from BUNIT if BUNIT is present
    hdu : int
        The index of the HDU to load into this spectrum.

    Notes
    -----
    Loader contributed by Kelle Cruz.
    """
    logging.info("Spectrum file looks like wcs1d-fits")

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
    spectrum: Spectrum1D
    file_name: str
        The path to the FITS file
    hdu: int
        Header Data Unit in FITS file to write to (base 0; default primary HDU)
    update_header: bool
        Update FITS header with all compatible entries in `spectrum.meta`
    unit : str or `~astropy.units.Unit`
        Unit for the flux (and associated uncertainty)
    dtype : str or `~numpy.dtype`
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


@data_loader('iraf', identifier=identify_iraf_wcs, dtype=Spectrum1D,
             extensions=['fits'])
def non_linear_wcs1d_fits(file_obj, spectral_axis_unit=None, flux_unit=None,
                          **kwargs):
    """Read wcs from files written by IRAF

    IRAF does not strictly follow the fits standard specially for non-linear
     wavelength solutions

    Parameters
    ----------

    file_obj: str, file-like or HDUList
        FITS file name, object (provided from name by Astropy I/O Registry),
        or HDUList (as resulting from astropy.io.fits.open()).

    spectral_axis_unit : `~astropy.Unit`, optional
        Spectral axis unit, default is None in which case will search for it
        in the header under the keyword 'WAT1_001'

    flux_unit : `~astropy.Unit`, optional
        Flux units, default is None. If not specified will attempt to read it
        using the keyword 'BUNIT' and if this keyword does not exist it will
        assume 'ADU'.

    Returns
    -------
    `specutils.Spectrum1D`
    """

    logging.info('Loading 1D non-linear fits solution')

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        header = hdulist[0].header
        for wcsdim in range(1, header['WCSDIM'] + 1):
            ctypen = header['CTYPE{:d}'.format(wcsdim)]
            if ctypen == 'LINEAR':
                logging.info("linear Solution: Try using "
                             "`format='wcs1d-fits'` instead")
                wcs = WCS(header)
                spectral_axis = _read_linear_iraf_wcs(wcs=wcs, dc_flag=header['DC-FLAG'])
            elif ctypen == 'MULTISPE':
                logging.info("Multi spectral or non-linear solution")
                spectral_axis = _read_non_linear_iraf_wcs(header=header, wcsdim=wcsdim)
            else:
                raise NotImplementedError

        if flux_unit is not None:
            data = hdulist[0].data * flux_unit
        elif 'BUNIT' in header:
            data = u.Quantity(hdulist[0].data, unit=header['BUNIT'])
        else:
            logging.info("Flux unit was not provided, neither it was in the"
                         "header. Assuming ADU.")
            data = u.Quantity(hdulist[0].data, unit='adu')

        if spectral_axis_unit is not None:
            spectral_axis *= spectral_axis_unit
        else:
            wat_head = header['WAT1_001']
            wat_dict = dict()
            for pair in wat_head.split(' '):
                wat_dict[pair.split('=')[0]] = pair.split('=')[1]
            if wat_dict['units'] == 'angstroms':
                logging.info("Found spectral axis units to be angstrom")
                spectral_axis *= u.angstrom

    meta = {'header': header}

    return Spectrum1D(flux=data, spectral_axis=spectral_axis, meta=meta)


def _read_linear_iraf_wcs(wcs, dc_flag):
    """Linear solution reader

    This method read the appropriate keywords. Calls the method _set_math_model
    which decides what is the appropriate mathematical model to be used and
    creates and then evaluates the model for an array.

    Parameters
    ----------

    wcs : `~astropy.wcs.WCS`
        Contains wcs information extracted from the header

    dc_flag : int
        Extracted from the header under the keyword DC-FLAG which defines what
        kind of solution is described. For linear solutions it is 0 or 1.

    Returns
    -------

    spectral_axis : `~numpy.ndarray`
        Mathematical model of wavelength solution evluated for each pixel
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

    header : `~astropy.io.fits.header.Header`
        Full header of file being loaded

    wcsdim : int
        Number of the wcs dimension to be read.

    Returns
    -------

    spectral_axis : `~numpy.ndarray`
        Mathematical model of wavelength solution evluated for each pixel
        position
    """

    wat_wcs_dict = {}
    ctypen = header['CTYPE{:d}'.format(wcsdim)]
    logging.info('Attempting to read CTYPE{:d}: {:s}'.format(wcsdim, ctypen))
    if ctypen == 'MULTISPE':
        # TODO (simon): What is the * (asterisc) doing here?.
        wat_head = header['WAT{:d}*'.format(wcsdim)]
        if len(wat_head) == 1:
            logging.debug('Get units')
            wat_array = wat_head[0].split(' ')
            for pair in wat_array:
                split_pair = pair.split('=')
                wat_wcs_dict[split_pair[0]] = split_pair[1]
                # print(wat_head[0].split(' '))
        elif len(wat_head) > 1:
            wat_string = ''
            for key in wat_head:
                wat_string += header[key]
            wat_array = shlex.split(wat_string.replace('=', ' '))
            if len(wat_array) % 2 == 0:
                for i in range(0, len(wat_array), 2):
                    # if wat_array[i] not in wcs_dict.keys():
                    wat_wcs_dict[wat_array[i]] = wat_array[i + 1]
                    # print(wat_array[i], wat_array[i + 1])

    for key in wat_wcs_dict.keys():
        logging.debug("{:d} -{:s}- {:s}".format(wcsdim,
                                                key,
                                                wat_wcs_dict[key]))

    if 'spec1' in wat_wcs_dict.keys():
        spec = wat_wcs_dict['spec1'].split()
        aperture = int(spec[0])
        beam = int(spec[1])
        disp_type = int(spec[2])
        disp_start = float(spec[3])
        disp_del_av = float(spec[4])
        pix_num = int(spec[5])
        dopp_fact = float(spec[6])
        aper_low = int(float(spec[7]))
        aper_high = int(float(spec[8]))
        weight = float(spec[9])
        zeropoint = float(spec[10])
        function_type = int(spec[11])
        order = int(float(spec[12]))
        min_pix_val = int(float(spec[13]))
        max_pix_val = int(float(spec[14]))

        params = [float(i) for i in spec[15:]]
        wcs_dict = {'aperture': aperture,
                    'beam': beam,
                    'dtype': disp_type,
                    'dstart': disp_start,
                    'avdelt': disp_del_av,
                    'pnum': pix_num,
                    'z': dopp_fact,
                    'alow': aper_low,
                    'ahigh': aper_high,
                    'weight': weight,
                    'zeropoint': zeropoint,
                    'ftype': function_type,
                    'order': order,
                    'pmin': min_pix_val,
                    'pmax': max_pix_val,
                    'fpar': params}

        logging.info('Retrieving model')
        math_model = _set_math_model(wcs_dict=wcs_dict)

        spectral_axis = math_model(range(1, wcs_dict['pnum'] + 1))
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
        return _none()
    elif wcs_dict['dtype'] == 0:
        return _linear_solution(wcs_dict=wcs_dict)
    elif wcs_dict['dtype'] == 1:
        return _log_linear(wcs_dict=wcs_dict)
    elif wcs_dict['dtype'] == 2:
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

        `~astropy.modeling.Model`

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
