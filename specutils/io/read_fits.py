# a simple fits readers that puts a lot of emphasis on the specwcs

from astropy.io import fits
from astropy import units as u


from specutils.wcs import specwcs
from specutils import Spectrum1D

def read_fits(filename, dispersion_unit=None, flux_unit=None):

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



