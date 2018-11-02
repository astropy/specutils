import astropy.units as u
from astropy.io import fits

from ...spectra import Spectrum1D
from ..registers import data_loader


def identify_jwst_fits(origin, *args, **kwargs):

    try:
        with fits.open(args[0]) as hdulist:
            # This is a near-guarantee that we have a JWST data product
            if not 'ASDF' in hdulist:
                return False
            # This indicates the data product contains  spectral data
            if not 'EXTRACT1D' in hdulist:
                return False
        return True
    # This probably means we didn't have a FITS file
    except Exception:
        return False


@data_loader("JWST", identifier=identify_jwst_fits, dtype=Spectrum1D)
def jwst_loader(filename, spectral_axis_unit=None, **kwargs):

    with fits.open(filename) as hdulist:
        # TODO: eventually we will pass a SpectrumCollection back, but for now
        # we just return the first spectrum in the file.
        hdu = hdulist['EXTRACT1D']
        wavelength = hdu.data['WAVELENGTH'] * u.Unit(hdu.header['TUNIT1'])
        flux = hdu.data['FLUX'] * u.Unit(hdu.header['TUNIT2'])
        error = hdu.data['ERROR'] * u.Unit(hdu.header['TUNIT3'])

        meta = dict(slitname=hdu.header.get('SLTNAME', ''))

    # TODO: pass uncertainty using the error from the HDU
    return Spectrum1D(flux=flux, spectral_axis=wavelength, meta=meta)
