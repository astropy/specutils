from astropy.io.fits import convenience, PrimaryHDU

multispec_cards = {'SIMPLE': True, 'BITPIX': -32, 'NAXIS': 2,
                   'CTYPE1': 'MULTISPE', 'CTYPE2': 'MULTISPE', 'WCSDIM': 2,
                   'EXTEND': False}

singlespec_cards = {'SIMPLE': True, 'BITPIX': -32, 'NAXIS': 1, 'WCSDIM': 1,
                   'EXTEND': False}

def write(spectrum, filename, clobber=True):
    if isinstance(spectrum, list):
        # assuming this list represents a mutispec system
        pass
    else:
        # assuming this spectrum has a polynomial WCS
        wcs = spectrum.wcs
        hdu = convenience._makehdu(spectrum.data, header=None)
        if hdu.is_image and not isinstance(hdu, PrimaryHDU):
            hdu = PrimaryHDU(spectrum.data, header=None)
        wcs.write_fits_header(hdu.header)

        hdu.writeto(filename, clobber=clobber)
