from astropy.io.fits import convenience, PrimaryHDU

multispec_cards = {'SIMPLE': True, 'BITPIX': -32, 'NAXIS': 2,
                   'CTYPE1': 'MULTISPE', 'CTYPE2': 'MULTISPE', 'WCSDIM': 2,
                   'EXTEND': False}

singlespec_cards = {'SIMPLE': True, 'BITPIX': -32, 'NAXIS': 1, 'WCSDIM': 1,
                   'EXTEND': False}

def write(filename, data, WCS, output_verify='exception', clobber=True,
          checksum=False):
    if len(data.shape) > 2:
        raise TypeError("data cannot have more than 2 dimensions")
    if len(data.shape) == 2 and \
            (not isinstance(WCS, list) or data.shape[1] != len(WCS)):
        raise TypeError("If data has more than one dimension then wcs must"
                        "be a list of the same length as the second dimension")

    hdu = convenience._makehdu(data, header=None)
    if hdu.is_image and not isinstance(hdu, PrimaryHDU):
        hdu = PrimaryHDU(data, header=None)

    if len(data.shape) == 2:
        hdu.header['CTYPE1'] = 'MULTISPE'
        hdu.header['CTYPE2'] = 'MULTISPE'
        hdu.header['WAT0_001'] = "system=multispec"
        unit_string = WCS[0].unit.to_string()
        if unit_string == "Angstrom":
            unit_string = "angstroms"
        label_string = "Wavelength"
        hdu.header['WAT1_001'] = "wtype=multispec label={0} unit={1}".\
            format(label_string, unit_string)
        specs = []
        for wcs in WCS:
            specs.append(wcs.get_fits_spec())
        # append spec to header
    else:
        WCS.write_fits_header(hdu.header)

    hdu.writeto(filename, clobber=clobber, output_verify=output_verify,
                checksum=checksum)
