from astropy.io.fits import convenience, PrimaryHDU
import numpy as np

multispec_cards = {'SIMPLE': True, 'BITPIX': -32, 'NAXIS': 2,
                   'CTYPE1': 'MULTISPE', 'CTYPE2': 'MULTISPE', 'WCSDIM': 2,
                   'EXTEND': False}

singlespec_cards = {'SIMPLE': True, 'BITPIX': -32, 'NAXIS': 1, 'WCSDIM': 1,
                   'EXTEND': False}


def _make_hdu(data, header=None):
    hdu = convenience._makehdu(data, header=header)
    if hdu.is_image and not isinstance(hdu, PrimaryHDU):
        hdu = PrimaryHDU(data, header=header)
    return hdu


def write(spectrum, filename, clobber=True):
    if isinstance(spectrum, list):
        # assuming this list represents a mutispec system, with each wcs an
        # IRAF combination WCS
        data = np.array([spectra.data for spectra in spectrum])
        hdu = _make_hdu(data)
        hdu.header['CTYPE1'] = 'MULTISPE'
        hdu.header['CTYPE2'] = 'MULTISPE'
        hdu.header['WCSDIM'] = 2
        hdu.header['WAT0_001'] = "system=multispec"
        unit_string = spectrum[0].wcs.unit.to_string()
        if unit_string == "Angstrom":
            unit_string = "angstroms"
        label_string = "Wavelength"
        hdu.header['WAT1_001'] = "wtype=multispec label={0} units={1}".\
                                    format(label_string, unit_string)
        spec_string = "wtype=multispec"
        for i, spectra in enumerate(spectrum):
            spec = " ".join(map(str, spectra.wcs.get_fits_spec()))
            spec_string += ' spec{0} = "{1}"'.format(i + 1, spec)
        wat_num = 1
        chars = 68
        for i in range(0, len(spec_string), chars):
            hdu.header["WAT2_{0:03d}".format(wat_num)] = spec_string[i: i+chars]
            wat_num += 1
    else:
        # assuming this spectrum has a polynomial WCS
        wcs = spectrum.wcs
        hdu = _make_hdu(spectrum.data)
        wcs.write_fits_header(hdu.header)

    hdu.writeto(filename, clobber=clobber)
