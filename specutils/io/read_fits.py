# a simple fits readers that puts a lot of emphasis on the specwcs
import re

from astropy.io import fits
from astropy import units as u
from astropy.utils import OrderedDict

from specutils.wcs import specwcs
from specutils import Spectrum1D

import numpy as np

wat_keyword_pattern = re.compile(
    '([^=\s]*)\s*=\s*(([^\"\'\s]+)|([\"\'][^\"\']+[\"\']))\s*')


class FITSWCSError(Exception):
    pass


class FITSWCSSpectrum1DError(FITSWCSError):
    pass


class FITSWCSSpectrum1DUnitError(FITSWCSError):
    pass


#keywords as described in http://iraf.net/irafdocs/specwcs.php

fits_wcs_spec_func_type = {1: 'chebyshev',
                           2: 'legendre',
                           3: 'cubicspline',
                           4: 'linearspline',
                           5: 'pixelcoordinatearray',
                           6: 'sampledcoordinatearray'}

wcs_attributes_function_parameters = {'chebyshev': ['order', 'pmin', 'pmax'],
                                      'legendre': ['order', 'pmin', 'pmax'],
                                      'linearspline': ['npieces', 'pmin',
                                                       'pmax'],
                                      'cubicspline': ['npieces', 'pmin', 'pmax']}

wcs_attributes_general_keywords = OrderedDict(
    [('aperture', int), ('beam', int), ('dispersion_type', int),
     ('dispersion0', float), ('average_dispersion_delta', float),
     ('no_valid_pixels', int), ('doppler_factor', float),
     ('aperture_low', float), ('aperture_high', float)])

wcs_attributes_function_keywords = OrderedDict(
    [('weight', float), ('zero_point_offset', float), ('type', int),
     ('order', int), ('pmin', float), ('pmax', float), ('npieces', int)])


def _get_num_coefficients(function_dict):
    """
    Returns the number of coeffecients according to the IRAF fits format defined
    here: http://iraf.net/irafdocs/specwcs.php
    """
    function_type = function_dict["type"]

    if function_type in ["legendre", "chebyshev"]:
        return function_dict["order"]
    elif function_type == "linearspline":
        return function_dict["npieces"] + 1
    elif function_type == "cubicspline":
        return function_dict["npieces"] + 3
    else:
        return 0


def _parse_fits_units(fits_unit_string):
    """
    Parse FITS units - only converting Angstroms to Angstrom

    Parameters
    ----------

    fits_unit_string: str
    """

    if fits_unit_string.lower().strip() == 'angstroms':
        fits_unit_string = 'Angstrom'

    return u.Unit(fits_unit_string)


def _parse_multispec_dict(multispec_dict):
    """
    Parse dictionary that contains the multispec information in wcs attributes (WAT keywords; often WAT2_???).
    each specN keywords has information as a string in the following format:
    ap beam dtype w1 dw nw z aplow aphigh wt_i w0_i ftype_i [parameters] [coefficients]

    Parameters
    ----------

        multispec_dict: dict-like object
            e.g. multi_spec_dict = {'wtype':'multispec', spec1:'...', 'spec2':'...', ..., 'specN':'...'}

    """

    parsed_multispec_dict = OrderedDict()

    for spec_key in multispec_dict:
        if not spec_key.lower().startswith('spec'):
            continue

        single_spec_dict = OrderedDict()
        spec_string = multispec_dict[spec_key].strip().split()
        for key_name, key_dtype in wcs_attributes_general_keywords.items():
            single_spec_dict[key_name] = key_dtype(spec_string.pop(0))
        single_spec_dict['functions'] = []
        while len(spec_string) > 0:
            # There seems to be a function defined for this spectrum -
            # checking that the dispersion type indicates that:
            assert single_spec_dict['dispersion_type'] == 2

            function_dict = {}
            for key_name, key_dtype in wcs_attributes_function_keywords.items():
                function_dict[key_name] = key_dtype(spec_string.pop(0))

                #last of the general keywords
                if key_name == 'type':
                    function_dict['type'] = fits_wcs_spec_func_type[function_dict['type']]
                    break

            # different function types defined in
            # http://iraf.net/irafdocs/specwcs.php -- see fits_wcs_spec_func_type

            function_type = function_dict['type']
            if function_type in wcs_attributes_function_parameters:
                for key_name in wcs_attributes_function_parameters[function_type]:
                    key_dtype = wcs_attributes_function_keywords[key_name]
                    function_dict[key_name] = key_dtype(spec_string.pop(0))
                num_coefficients = _get_num_coefficients(function_dict)
                coefficients = spec_string[:num_coefficients]
                function_dict['coefficients'] = list(map(float, coefficients))
                spec_string = spec_string[num_coefficients:]
            else:
                raise NotImplementedError

            single_spec_dict["functions"].append(function_dict)

        parsed_multispec_dict[spec_key] = single_spec_dict

    return parsed_multispec_dict


class FITSWCSSpectrum(object):
    """This class is designed to extract all known spectroscopic WCS keywords from a FITS keyword header.
     The resulting keywords can then be validated (several keywords will encode the same information). Other keywords
     come in pairs (e.g. CRVAL1, CRVAL2) and will be available as a list (e.g. `.affine_transform_dict['crval']`).
     Several FITS readers make use of this object.

    Parameters
    ----------

    fits_header: dict-like object (e.g. `~astropy.io.fits.Header`)
        FITS Header to be read

    """

    def __init__(self, fits_header):
        self.fits_header = fits.Header(fits_header)

        self.naxis = self.fits_header['naxis']

        self.shape = []

        for i in range(self.naxis):
            self.shape.append(self.fits_header['naxis{0:d}'.format(i + 1)])

        try:
            self.wcs_dim = self.fits_header['WCSDIM']
        except KeyError:
            self.wcs_dim = None

        try:
            self.global_wcs_attributes = self.read_wcs_attributes(0)
        except FITSWCSError:
            self.global_wcs_attributes = None

        if self.wcs_dim is not None:
            self.wcs_attributes = []
            for axis in range(self.wcs_dim):
                self.wcs_attributes.append(self.read_wcs_attributes(axis + 1))

        self.affine_transform_dict, self.transform_matrix = self.read_affine_transforms()
        self.units = self.read_wcs_units()


    def read_affine_transforms(self, wcs_dim=None):

        if wcs_dim is None:
            if self.wcs_dim is None:
                wcs_dim = self.fits_header['NAXIS']
            else:
                wcs_dim = self.wcs_dim

        affine_transform_keywords = ('ctype', 'crpix', 'crval', 'cdelt')
        affine_transform_dict = dict(
            [(key, [None] * wcs_dim) for key in affine_transform_keywords])

        for i in range(wcs_dim):
            for key in affine_transform_dict:
                affine_transform_dict[key][i] = self.fits_header.get(
                    '{0:s}{1:d}'.format(key, i + 1))

        transform_matrix = self.read_transform_matrix(wcs_dim, cdelt=
        affine_transform_dict['cdelt'])

        return affine_transform_dict, transform_matrix

    def read_wcs_units(self, wcs_dim=None):

        if wcs_dim is None:
            if self.wcs_dim is None:
                wcs_dim = self.fits_header['NAXIS']
            else:
                wcs_dim = self.wcs_dim

        units = [None] * wcs_dim

        for i in range(wcs_dim):
            try:
                cunit_string = self.fits_header['cunit{0:d}'.format(i + 1)]
            except KeyError:
                continue

            if cunit_string.strip().lower() == 'angstroms':
                units[i] = u.AA
            else:
                units[i] = u.Unit(cunit_string)

        return units


    def read_transform_matrix(self, matrix_dim, cdelt=None):
        if len(self.fits_header['cd?_?']) > 0:
            if matrix_dim is None:
                if self.wcs_dim is None:
                    matrix_dim = self.fits_header['NAXIS']
                else:
                    matrix_dim = self.wcs_dim

            transform_matrix = np.matrix(np.zeros((matrix_dim, matrix_dim)))
            matrix_element_keyword_pattern = re.compile('cd(\d)_(\d)',
                                                        re.IGNORECASE)
            for matrix_element_keyword in self.fits_header['cd?_?']:
                i, j = map(int, matrix_element_keyword_pattern.match(
                    matrix_element_keyword).groups())
                transform_matrix[j - 1, i - 1] = self.fits_header[
                    matrix_element_keyword]
            if cdelt is not None:
                for i in range(matrix_dim):
                    if cdelt[i] is not None:
                        np.testing.assert_almost_equal(cdelt[i],
                                                       transform_matrix[i, i])

            return transform_matrix

        else:
            return None


    def read_wcs_attributes(self, axis):
        """
        Reading WCS attribute information in WAT0_001-like keywords

        Parameters
        ----------

        axis: int
            specifying which axis to read (e.g axis=2 will read WAT2_???).
        """

        wcs_attributes = self.fits_header['wat{0:d}_???'.format(axis)]
        if len(wcs_attributes) == 0:
            raise FITSWCSError

        raw_wcs_attributes = ''.join([wcs_attributes[key].ljust(68) for key in
                                      sorted(wcs_attributes.keys())])

        wat_dictionary = OrderedDict()
        for wat_keyword_match in wat_keyword_pattern.finditer(
                raw_wcs_attributes):
            wat_dictionary[wat_keyword_match.groups()[0]] = \
                wat_keyword_match.groups()[1].strip('\"\'')

        if 'units' in wat_dictionary:
            wat_dictionary['units'] = _parse_fits_units(wat_dictionary['units'])

        return wat_dictionary

def multispec_wcs_reader(wcs_info, dispersion_unit=None):
    """Extracting multispec information out of WAT header keywords and
    building WCS with it

    Parameters
    ----------

    dispersion_unit : astropy.unit.Unit, optional
        specify a unit for the dispersion if none exists or overwrite,
        default=None
    """

    assert wcs_info.global_wcs_attributes['system'] == 'multispec'
    assert wcs_info.wcs_attributes[1]['wtype'] == 'multispec'

    if dispersion_unit is None:
        dispersion_unit = wcs_info.wcs_attributes[0]['units']

    multispec_dict = _parse_multispec_dict(wcs_info.wcs_attributes[1])
    wcs_dict = OrderedDict()
    for spec_key in multispec_dict:
        single_spec_dict = multispec_dict[spec_key]
        if single_spec_dict['dispersion_type'] == 1:
            #log-linear dispersion
            log = True
        else:
            log = False

        if single_spec_dict['dispersion_type'] in [0, 1]:
            #linear or log-linear dispersion
            dispersion_wcs = specwcs.Spectrum1DPolynomialWCS(
                degree=1, c0=single_spec_dict["dispersion0"],
                c1=single_spec_dict["average_dispersion_delta"])

        else:
            # single_spec_dict['dispersion_type'] == 2
            dispersion_wcs = specwcs.WeightedCombinationWCS()
            for function_dict in single_spec_dict["functions"]:
                if function_dict['type'] == 'legendre':
                    ##### @embray can you figure out if that's the only way to
                    ##  instantiate a polynomial (with c0=xx, c1=xx, ...)?

                    coefficients = dict([('c{:d}'.format(i),
                                        function_dict['coefficients'][i])
                                        for i in range(function_dict['order'])])

                    wcs = specwcs.Spectrum1DIRAFLegendreWCS(
                        function_dict['order'], function_dict['pmin'],
                        function_dict['pmax'], **coefficients)

                elif function_dict['type'] == 'chebyshev':
                    coefficients = dict([('c{:d}'.format(i),
                                        function_dict['coefficients'][i])
                                        for i in range(function_dict['order'])])

                    wcs = specwcs.Spectrum1DIRAFChebyshevWCS(
                        function_dict['order'], function_dict['pmin'],
                        function_dict['pmax'], **coefficients)

                elif function_dict['type'] in ['linearspline', 'cubicspline']:
                    if function_dict['type'] == 'linearspline':
                        degree = 1
                    else:
                        degree = 3
                    n_pieces = function_dict['npieces']
                    pmin = function_dict['pmin']
                    pmax = function_dict['pmax']
                    y = [function_dict['coefficients'][i]
                         for i in range(n_pieces + degree)]
                    wcs = specwcs.Spectrum1DIRAFBSplineWCS(degree, n_pieces, y,
                                                           pmin, pmax)
                else:
                    raise NotImplementedError
                dispersion_wcs.add_WCS(
                    wcs, weight=function_dict["weight"],
                    zero_point_offset=function_dict["zero_point_offset"])

        composite_wcs = specwcs.MultispecIRAFCompositeWCS(
            dispersion_wcs, single_spec_dict["no_valid_pixels"],
            z=single_spec_dict["doppler_factor"], log=log,
            aperture=single_spec_dict["aperture"], beam=single_spec_dict["beam"]
            , aperture_low=single_spec_dict["aperture_low"],
            aperture_high=single_spec_dict["aperture_high"],
            unit=dispersion_unit)
        wcs_dict[spec_key] = composite_wcs
    return wcs_dict


def read_fits_wcs_linear1d(wcs_info, dispersion_unit=None, spectral_axis=0):
    """Read very a very simple 1D WCS mainly comprising of CRVAL, CRPIX, ...
    from a FITS WCS Information container

    Parameters
    ----------

    fits_wcs_information : ~specutils.io.read_fits.FITSWCSSpectrum
        object compiling WCS information to be used in these readers

    .. TODO: Examine this further.

    """

    # for the 1D reader setting the spectral_axis to anything else than 0 seems
    # to be strange;
    # actually, it's perfectly reasonable IF you want to extract
    # information from, say, a data cube and you're only interested in the
    # spectral dimension.  This tool can easily be used for that purpose.

    dispersion_unit = dispersion_unit

    dispersion_delta = None
    if wcs_info.transform_matrix is not None:
        dispersion_delta = wcs_info.transform_matrix[spectral_axis, spectral_axis]

    if wcs_info.affine_transform_dict['cdelt'][spectral_axis] is not None:
        #checking that both cd1_1 and cdelt1 are either the same or
        # one of them non-existent
        if dispersion_delta is not None:
            np.testing.assert_almost_equal(
                dispersion_delta,
                wcs_info.affine_transform_dict['cdelt'][spectral_axis])
        dispersion_delta = wcs_info.affine_transform_dict['cdelt'][spectral_axis]

    if dispersion_delta is None:
        raise FITSWCSSpectrum1DError

    if wcs_info.affine_transform_dict['crval'][spectral_axis] is None:
        raise FITSWCSSpectrum1DError
    else:
        dispersion_start = wcs_info.affine_transform_dict['crval'][spectral_axis]

    pixel_offset = wcs_info.affine_transform_dict['crpix'][spectral_axis] or 1
    pixel_offset -= 1

    dispersion_unit = wcs_info.units[spectral_axis] or dispersion_unit

    if None in [dispersion_start, dispersion_delta, pixel_offset]:
        raise FITSWCSSpectrum1DError
    dispersion_start += -pixel_offset * dispersion_delta
    return specwcs.Spectrum1DPolynomialWCS(degree=1, unit=dispersion_unit,
                                           c0=dispersion_start,
                                           c1=dispersion_delta)

def read_fits_spectrum1d(filename, dispersion_unit=None, flux_unit=None):
    """
    1D reader for spectra in FITS format. This function determines what format
    the FITS file is in, and attempts to read the Spectrum. This reader just
    uses the primary extension in a FITS file and reads the data and header from
    that. It will return a Spectrum1D object if the data is linear, or a list of
    Spectrum1D objects if the data format is multi-spec

    Parameters
    ----------

    filename : str
        FITS filename

    dispersion_unit : ~astropy.unit.Unit, optional
        unit of the dispersion axis - will overwrite possible information given
        in the FITS keywords
        default = None

    flux_unit : ~astropy.unit.Unit, optional
        unit of the flux

    Raises
    --------
    NotImplementedError
        If the format can't be read currently
    """
    if dispersion_unit:
        dispersion_unit = u.Unit(dispersion_unit)

    data = fits.getdata(filename)
    header = fits.getheader(filename)

    wcs_info = FITSWCSSpectrum(header)

    if wcs_info.naxis == 1:
        wcs = read_fits_wcs_linear1d(wcs_info, dispersion_unit=dispersion_unit)
        return Spectrum1D(data, wcs=wcs, unit=flux_unit)
    elif wcs_info.naxis == 2 and \
            wcs_info.affine_transform_dict['ctype'] == ["MULTISPE", "MULTISPE"]:
        multi_wcs = multispec_wcs_reader(wcs_info, dispersion_unit=dispersion_unit)
        multispec = []
        for spectrum_data, spectrum_wcs in zip(data, multi_wcs.values()):
            multispec.append(
                Spectrum1D(spectrum_data, wcs=spectrum_wcs, unit=flux_unit))
        return multispec

    elif wcs_info.naxis == 3 and \
            wcs_info.affine_transform_dict['ctype'] == ["LINEAR","LINEAR","LINEAR"]:
        wcs = read_fits_wcs_linear1d(wcs_info, dispersion_unit=dispersion_unit)
        equispec = []
        for i in range(data.shape[0]):
            equispec.append(
                Spectrum1D(data[i][0], wcs=wcs, unit=flux_unit))
        return equispec
        
    elif wcs_info.naxis == 3 and \
            wcs_info.affine_transform_dict['ctype'] == ["MULTISPE", "MULTISPE","LINEAR"]:
        multi_wcs = multispec_wcs_reader(wcs_info, dispersion_unit=dispersion_unit)
        multispec = []
        for j in range(data.shape[1]):
            equispec = []
            for i in range(data.shape[0]):
                equispec.append(
                    Spectrum1D(data[i][j], wcs=list(multi_wcs.values())[j], unit=flux_unit))
            multispec.append(equispec)
        return multispec

    else:
        raise NotImplementedError("Either the FITS file does not represent a 1D"
                                  " spectrum or the format isn't supported yet")

