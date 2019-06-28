import astropy.units as u
from astropy.io import fits

from ...spectra import Spectrum1D, SpectrumList
from ..registers import data_loader


def identify_jwst_fits(origin, *args, **kwargs):
    """
    Check whether the given file is a JWST spectral data product.

    This check is fairly simple. It expects FITS files that contain an ASDF
    header (which is not used here, but indicates a JWST data product). It then
    looks for at least one EXTRACT1D header, which contains spectral data.
    """

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


@data_loader("JWST", identifier=identify_jwst_fits, dtype=SpectrumList,
             extensions=['fits'])
def jwst_loader(filename, spectral_axis_unit=None, **kwargs):
    """
    Loader for JWST data files.

    Parameters
    ----------
    file_name: str
        The path to the FITS file

    Returns
    -------
    data: SpectrumList
        A list of the spectra that are contained in this file.
    """

    spectra = []

    with fits.open(filename) as hdulist:
        for hdu in hdulist:
            if hdu.name != 'EXTRACT1D':
                continue

            # Provide reasonable defaults based on the units assigned by the
            # extract1d step of the JWST pipeline. TUNIT fields should be
            # populated by the pipeline, but it's apparently possible for them
            # to be missing in some files.
            wavelength_units = u.Unit(hdu.header.get('TUNIT1', 'um'))
            flux_units = u.Unit(hdu.header.get('TUNIT2', 'mJy'))
            error_units = u.Unit(hdu.header.get('TUNIT3', 'mJy'))

            wavelength = hdu.data['WAVELENGTH'] * wavelength_units
            flux = hdu.data['FLUX'] * flux_units
            error = hdu.data['ERROR'] * error_units

            meta = dict(slitname=hdu.header.get('SLTNAME', ''))

            # TODO: pass uncertainty using the error from the HDU
            spec = Spectrum1D(flux=flux, spectral_axis=wavelength, meta=meta)
            spectra.append(spec)

    return SpectrumList(spectra)
