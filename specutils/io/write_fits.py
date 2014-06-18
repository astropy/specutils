from astropy.io.fits import convenience, PrimaryHDU

multispec_cards = {'SIMPLE': True, 'BITPIX': -32, 'NAXIS': 2,
                   'CTYPE1': 'MULTISPE', 'CTYPE2': 'MULTISPE', 'WCSDIM': 2,
                   'EXTEND': False}

singlespec_cards = {'SIMPLE': True, 'BITPIX': -32, 'NAXIS': 1, 'WCSDIM': 1,
                   'EXTEND': False}

def write(filename, data, wcs, output_verify='exception', clobber=True,
          checksum=False):
    if len(data.shape) > 2:
        raise TypeError("data cannot have more than 2 dimensions")
    if len(data.shape) == 2 and \
            (not isinstance(wcs, list) or data.shape[1] != len(wcs)):
        raise TypeError("If data has more than one dimension then wcs must"
                        "be a list of the same length as the second dimension")

    hdu = convenience._makehdu(data, header=None)
    if hdu.is_image and not isinstance(hdu, PrimaryHDU):
        hdu = PrimaryHDU(data, header=None)

    if len(data.shape) == 2:
        for wcs1D in wcs:
            wcs1D.add_to_header(hdu.header)
    else:
        wcs.write_fits_header(hdu.header)

    hdu.writeto(filename, clobber=clobber, output_verify=output_verify,
                checksum=checksum)
