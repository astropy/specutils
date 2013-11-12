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


class BaseReadFITSSpectrumWCS(object):
    """

    """

    def __init__(self, fits_header):
        self.fits_header = fits.Header(fits_header)
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


            if self.wcs_dim is None:
                matrix_dim = self.fits_header['NAXIS']
            else:
                matrix_dim = self.wcs_dim

        self.transform_matrix = self.read_transform_matrix()



    def read_affine_transforms(self, wcs_dim=None):
        if self.wcs_dim is none


    def read_transform_matrix(self, matrix_dim=None):
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

        raw_wcs_attributes = ''.join([wcs_attributes[key] for key in sorted(wcs_attributes.keys())])

        wat_dictionary = OrderedDict()
        for wat_keyword_match in wat_keyword_pattern.finditer(raw_wcs_attributes):
            wat_dictionary[wat_keyword_match.groups()[0]] = wat_keyword_match.groups()[1].strip('\"\'')

        return wat_dictionary






class LinearReadFITSSpectrumWCS(BaseReadFITSSpectrumWCS):

    def _init__(self, fits_header):
        super(LinearReadFITSSpectrumWCS, self).__init__(fits_header)






def read_fits_single_spec(filename, dispersion_unit=None, flux_unit=None):

    if dispersion_unit:
        dispersion_unit = u.Unit(dispersion_unit)
    data = fits.getdata(filename)
    header = fits.getheader(filename)

    for fits_wcs in specwcs.fits_capable_wcs:
        try:
            wcs = fits_wcs.from_fits_header(header, unit=dispersion_unit)
        except specwcs.Spectrum1DWCSFITSError:
            continue
        except specwcs.Spectrum1DWCSUnitError:
            raise specwcs.Spectrum1DWCSUnitError('%s can read WCS information in the file, however no dispersion unit'
                                                 ' was found. Please specify this using the keyword dispersion_unit '
                                                 '(e.g. dispersion_unit=\'Angstrom\')' % fits_wcs)
        else:
            break
    else:
        raise specwcs.Spectrum1DWCSError('File %s does not contain a specutils-readable WCS - exiting' % filename)

    return Spectrum1D(data, wcs=wcs, unit=flux_unit)


def combine_headers(headers):
    combined = []
    for header in headers.values():
        # be sure to append back any trailing spaces
        combined.append('{:68s}'.format(header))
    combined = "".join(combined)
    combined = combined.replace('= ', '=')
    combined = combined.replace(' =', '=')
    return combined


def read_fits_multispec(filename, dispersion_unit=None, flux_unit=None):

    if dispersion_unit:
        dispersion_unit = u.Unit(dispersion_unit)

    data = fits.getdata(filename)
    header = fits.getheader(filename)

    try:
        assert header['WAT0_001'] == 'system=multispec'
    except(KeyError, AssertionError):
        raise specwcs.Spectrum1DWCSFITSError

    head1 = combine_headers(header['WAT1_*'])
    head1_unit = re.match(r'.*units+=(.*)', head1).groups()[0].strip()
    multispec_unit = u.Unit(head1_unit.lower().replace('angstroms',
                                                       'Angstrom'))
    if dispersion_unit:
        scale = multispec_unit.to(dispersion_unit)
    else:
        scale = 1
        dispersion_unit = multispec_unit

    head2 = combine_headers(header['WAT2_*'])
    specs = re.findall(r'spec[0-9]+="[0-9eE .+-]*"', head2)
    multispec = []

    for data1, spec in zip(data, specs):
        spec = spec[4:].replace('=',' ').replace('"','')
        #specN = ap beam dtype w1 dw nw z aplow aphigh [functions_i]
        #    0    1  2     3    4  5  6 7   8     9     10:
        spec = spec.split()
        specwcs_keys = {}
        for key, fmt in (('N', int), ('ap', int), ('beam', int),
                         ('dtype', int), ('w1', float), ('dw', float),
                         ('nw', int), ('z', float), ('aplow', float),
                         ('aphigh', float)):
            specwcs_keys[key] = fmt(spec.pop(0))

        if specwcs_keys['dtype'] == -1:   # no dispersion solution in header
            pass

        elif specwcs_keys['dtype'] == 0:  # linear dispersion
            # TBD: do something with aplow, aphigh
            wcs = specwcs_keys.Spectrum1DLinearWCS(specwcs_keys['w1'],
                                                   specwcs_keys['dw'], 0)

        elif specwcs_keys['dtype'] == 1:  # log-linear dispersion
            raise NotImplementedError

        elif specwcs_keys['dtype'] == 2:  # non-linear dispersion
            while len(spec) > 0:
                function = {}
                function['weight'] = float(spec.pop(0))
                function['offset'] = float(spec.pop(0))
                func = int(spec.pop(0))
                if func < 5:
                    if func == 1:  # Chebyshev
                        return NotImplementedError

                    elif func == 2:
                        npol = int(spec.pop(0))
                        apmin = float(spec.pop(0))
                        apmax = float(spec.pop(0))
                        coeff = dict([('c{:d}'.format(i),
                                       float(spec.pop(0)) * scale)
                                      for i in range(npol)])
                        wcs = specwcs.Spectrum1DLegendreWCS(
                            npol - 1, domain=[apmin, apmax],
                            unit=dispersion_unit, **coeff)

                    elif func == 3:  # linear spline
                        raise NotImplementedError

                    elif func == 4:  # cubic spline
                        raise NotImplementedError

                    elif func == 5:  # pixel coordinate array
                        raise NotImplementedError

                    elif func == 6:  # sampled coordinate array
                        raise NotImplementedError
                    else:
                        raise specwcs.Spectrum1DWCSFITSError(
                            'File does not conform to IRAF format - '
                            'non-linear of dispersion function must be '
                            '[1,2,3,4,5,6]')
                if len(spec):
                    raise specwcs.Spectrum1DWCSFITSError(
                        'Cannot yet deal with multiple functions')

        else:
            raise specwcs.Spectrum1DWCSFITSError(
                'File does not conform to IRAF format - '
                'type of dispersion function must be [-1,0,1,2]')

        multispec += [Spectrum1D(data1, wcs=wcs)]

    return multispec
