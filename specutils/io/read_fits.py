# a simple fits readers that puts a lot of emphasis on the specwcs
import re

from astropy.io import fits
from astropy import units as u
from collections import OrderedDict

from specutils.wcs import specwcs
from specutils import Spectrum1D

import numpy as np

wat_keyword_pattern = re.compile('([^=\s]*)\s*=\s*(([^\"\'\s]+)|([\"\'][^\"\']+[\"\']))\s*')

class FITSWCSError(Exception):
    pass

class FITSWCSSpectrum1DError(FITSWCSError):
    pass

class FITSWCSSpectrum1DUnitError(FITSWCSError):
    pass



fits_wcs_spec_func_type = {1:'chebyshev',
                           2:'legendre',
                           3:'cubicspline',
                           4:'linearspline',
                           5:'pixelcoordinatearray',
                           6:'sampledcoordinatearray'
        }

fits_wcs_spec_func_parameter = {'chebyshev': ['order', 'pmin', 'pmax'],
                                'legendre' : ['order', 'pmin', 'pmax']}

fits_wcs_spec_func_parameter_dtype = {'chebyshev': [int, float, float],
                                'legendre' : [int, float, float]}

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
    general_spec_key_name = ('ap', 'beam', 'dtype', 'w1', 'dw', 'nw', 'z', 'aplow', 'aphigh')
    general_spec_key_dtype = (int, int, int, float, float, float, float, float, float)

    function_spec_key_name = ('weight', 'zero_point_offset', 'type')
    function_spec_key_dtype = (float, float, int)

    parsed_multispec_dict = OrderedDict()

    for spec_key in multispec_dict:
        if not spec_key.lower().startswith('spec'):
            continue

        single_spec_dict = OrderedDict()
        split_single_spec_string = multispec_dict[spec_key].strip().split()
        for key_name, key_dtype in zip(general_spec_key_name, general_spec_key_dtype):
            single_spec_dict[key_name] = key_dtype(split_single_spec_string.pop(0))

        if len(split_single_spec_string) > 0:

            #There seems to be a function defined for this spectrum - checking that the dispersion type indicates that:
            assert single_spec_dict['dtype'] == 2

            single_spec_dict['function'] = {}
            for key_name, key_dtype in zip(function_spec_key_name, function_spec_key_dtype):
                single_spec_dict['function'][key_name] = key_dtype(split_single_spec_string.pop(0))
                if key_name == 'type':
                    single_spec_dict['function']['type'] = fits_wcs_spec_func_type[single_spec_dict['function']['type']]
            #different function types defined in http://iraf.net/irafdocs/specwcs.php -- see fits_wcs_spec_func_type

            function_type =  single_spec_dict['function']['type']
            if function_type in fits_wcs_spec_func_parameter:

                for key_name, key_dtype in zip(fits_wcs_spec_func_parameter[function_type],
                                               fits_wcs_spec_func_parameter_dtype[function_type]):
                    single_spec_dict['function'][key_name] = key_dtype(split_single_spec_string.pop(0))

                single_spec_dict['function']['coefficients'] = map(float, split_single_spec_string)

            else:
                raise NotImplementedError

        parsed_multispec_dict[spec_key] = single_spec_dict

    return parsed_multispec_dict








class FITSWCSSpectrum(object):
    """

    """

    def __init__(self, fits_header):
        self.fits_header = fits.Header(fits_header)

        self.naxis = self.fits_header['naxis']

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
            for axis in xrange(self.wcs_dim):
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
        affine_transform_dict = dict([(key, [None] * wcs_dim) for key in affine_transform_keywords])

        for i in xrange(wcs_dim):
            for key in affine_transform_dict:
                affine_transform_dict[key][i] = self.fits_header.get('{0:s}{1:d}'.format(key, i+1))

        transform_matrix = self.read_transform_matrix(wcs_dim, cdelt=affine_transform_dict['cdelt'])

        return affine_transform_dict, transform_matrix

    def read_wcs_units(self, wcs_dim=None):

        if wcs_dim is None:
            if self.wcs_dim is None:
                wcs_dim = self.fits_header['NAXIS']
            else:
                wcs_dim = self.wcs_dim

        units = [None] * wcs_dim

        for i in xrange(wcs_dim):
            try:
                cunit_string = self.fits_header['cunit{0:d}'.format(i+1)]
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
            matrix_element_keyword_pattern = re.compile('cd(\d)_(\d)', re.IGNORECASE)
            for matrix_element_keyword in self.fits_header['cd?_?']:
                i, j = map(int, matrix_element_keyword_pattern.match(matrix_element_keyword).groups())
                transform_matrix[j-1, i-1] = self.fits_header[matrix_element_keyword]
            if cdelt is not None:
                for i in xrange(matrix_dim):
                    if cdelt[i] is not None:
                        np.testing.assert_almost_equal(cdelt[i], transform_matrix[i,i])

            return transform_matrix

        else:
            return None


    def read_wcs_attributes(self, axis):
        """
        Reading WCS attribute information in WAT0_001-like keywords

        Paramaters:
            axis: int
                specifying which axis to read (e.g axis=2 will read WAT2_???).
        """
        wcs_attributes = self.fits_header['wat{0:d}_???'.format(axis)]
        if len(wcs_attributes) == 0:
            raise FITSWCSError

        raw_wcs_attributes = ''.join([wcs_attributes[key].ljust(68) for key in sorted(wcs_attributes.keys())])

        wat_dictionary = OrderedDict()
        for wat_keyword_match in wat_keyword_pattern.finditer(raw_wcs_attributes):
            wat_dictionary[wat_keyword_match.groups()[0]] = wat_keyword_match.groups()[1].strip('\"\'')

        if 'units' in wat_dictionary:
            wat_dictionary['units'] = _parse_fits_units(wat_dictionary['units'])

        return wat_dictionary

    def get_multispec_wcs(self, dispersion_unit=None):
        """
            get multispec_wcs out of the existing information
        """

        assert self.naxis == 2
        assert self.global_wcs_attributes['system'] == 'multispec'
        assert self.wcs_attributes[1]['wtype'] == 'multispec'

        if dispersion_unit is None:
            dispersion_unit = self.wcs_attributes[0]['units']

        multispec_dict = _parse_multispec_dict(self.wcs_attributes[1])
        multispec_wcs_dict = OrderedDict()
        for spec_key in multispec_dict:
            single_spec_dict = multispec_dict[spec_key]

            if single_spec_dict['function']['type'] == 'legendre':
                function_dict = single_spec_dict['function']

                ##### @embray can you figure out if that's the only way to instantiate a polynomial (with c0=xx, c1=xx, ...)?

                coefficients = dict([('c{:d}'.format(i), function_dict['coefficients'][i])
                                     for i in range(function_dict['order'])])
                multispec_wcs_dict[spec_key] = specwcs.Spectrum1DLegendreWCS(function_dict['order'] - 1,
                                                                   domain=[function_dict['pmin'],
                                                                           function_dict['pmax']],
                                                                   unit=dispersion_unit,
                                                                   **coefficients)


            else:
                raise NotImplementedError
        return multispec_wcs_dict






def read_fits_wcs_linear1d(fits_wcs_information, dispersion_unit=None):

    dispersion_unit = dispersion_unit

    if fits_wcs_information.naxis != 1:
        raise FITSWCSSpectrum1DError

    dispersion_delta = None
    if fits_wcs_information.transform_matrix is not None:
        dispersion_delta = fits_wcs_information.transform_matrix[0, 0]


    if fits_wcs_information.affine_transform_dict['cdelt'][0] is not None:
        #checking that both cd1_1 and cdelt1 are either the same or one of them non-existent
        if dispersion_delta is not None:
            assert np.testing.assert_almost_equal(dispersion_delta, fits_wcs_information.affine_transform_dict['cdelt1'])
        dispersion_delta = fits_wcs_information.affine_transform_dict['cdelt'][0]

    if dispersion_delta is None:
        raise FITSWCSSpectrum1DError

    if fits_wcs_information.affine_transform_dict['crval'][0] is None:
        raise FITSWCSSpectrum1DError
    else:
        dispersion_start = fits_wcs_information.affine_transform_dict['crval'][0]

    pixel_offset = fits_wcs_information.affine_transform_dict['crpix'][0] or 1
    pixel_offset -= 1


    dispersion_unit = fits_wcs_information.units[0] or dispersion_unit


    if None in [dispersion_start, dispersion_delta, pixel_offset]:
        raise FITSWCSSpectrum1DError

    if dispersion_unit is None:
        raise FITSWCSSpectrum1DUnitError
    return specwcs.Spectrum1DLinearWCS(dispersion_start, dispersion_delta, pixel_offset, dispersion_unit)



spectrum1d_wcs_readers = [read_fits_wcs_linear1d]


def read_fits_spectrum1d(filename, dispersion_unit=None, flux_unit=None):

    if dispersion_unit:
        dispersion_unit = u.Unit(dispersion_unit)

    data = fits.getdata(filename)
    header = fits.getheader(filename)

    fits_wcs_information = FITSWCSSpectrum(header)

    for fits_wcs in spectrum1d_wcs_readers:
        try:
            wcs = fits_wcs(fits_wcs_information, dispersion_unit=dispersion_unit)
        except FITSWCSSpectrum1DError:
            continue
        except FITSWCSSpectrum1DUnitError:
            raise specwcs.Spectrum1DWCSUnitError('"{0:s}" can read WCS information in the file, however no dispersion unit'
                                                 ' was found. Please specify this using the keyword dispersion_unit '
                                                 '(e.g. dispersion_unit=\'Angstrom\')'.format(fits_wcs.__name__))
        else:
            break
    else:
        raise specwcs.Spectrum1DWCSError('File {0:s} does not contain a spectrum1d readable WCS'.format(filename))

    return Spectrum1D(data, wcs=wcs, unit=flux_unit)



def read_fits_multispec_to_list(filename, dispersion_unit=None, flux_unit=None):

    if dispersion_unit:
        dispersion_unit = u.Unit(dispersion_unit)



    data = fits.getdata(filename)
    header = fits.getheader(filename)

    fits_wcs_information = FITSWCSSpectrum(header)

    multispec_wcs = fits_wcs_information.get_multispec_wcs(dispersion_unit)

    multispec = []
    for spectrum_data, spectrum_wcs in zip(data, multispec_wcs.values()):
        multispec.append(Spectrum1D(spectrum_data, wcs=spectrum_wcs))

    return multispec
