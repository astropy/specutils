import warnings

import astropy.units as u
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.utils.exceptions import AstropyUserWarning
import numpy as np

from ...spectra import Spectrum1D, SpectrumList
from ..registers import data_loader


__all__ = ["jwst_loader"]


def identify_jwst_x1d_fits(origin, *args, **kwargs):
    """
    Check whether the given file is a JWST x1d spectral data product.
    """

    try:
        with fits.open(args[0]) as hdulist:
            # This is a near-guarantee that we have a JWST data product
            if not 'ASDF' in hdulist:
                return False
            # This indicates the data product contains  spectral data
            if not 'EXTRACT1D' in hdulist:
                return False
            if not hdulist[0].header["TELESCOP"] == "JWST":
                return False
        return True
    # This probably means we didn't have a FITS file
    except Exception:
        return False


@data_loader("JWST", identifier=identify_jwst_x1d_fits, dtype=SpectrumList,
             extensions=['fits'])
def jwst_loader(filename, **kwargs):
    """
    Loader for JWST x1d 1-D spectral data in FITS format

    Parameters
    ----------
    filename: str
        The path to the FITS file

    Returns
    -------
    SpectrumList
        A list of the spectra that are contained in the file.
    """

    spectra = []

    with fits.open(filename) as hdulist:

        for hdu in hdulist:
            # Read only the BinaryTableHDUs named EXTRACT1D
            if hdu.name != 'EXTRACT1D':
                continue

            wavelength_units = u.Unit(hdu.columns["wavelength"].unit)
            wavelength = hdu.data["wavelength"] * wavelength_units

            # Determine if FLUX or SURF_BRIGHT column should be returned
            # based on whether it is point or extended source
            try:
                srctype = hdu.header.get("srctype")
            except KeyError:
                # SRCTYPE is in primary header because there is only one spectrum
                srctype = hdulist.header.get("srctype")

            if srctype == "POINT":
                flux_units = u.Unit(hdu.columns["flux"].unit)
                flux = hdu.data["flux"] * flux_units

                error_units = u.Unit(hdu.columns["error"].unit)
                uncertainty = StdDevUncertainty(hdu.data["error"] * error_units)

            elif srctype == "EXTENDED":
                flux_units = u.Unit(hdu.columns["surf_bright"].unit)
                flux = hdu.data["surf_bright"] * flux_units

                error_units = u.Unit(hdu.columns["sb_error"].unit)
                uncertainty = StdDevUncertainty(hdu.data["sb_error"] * error_units)

            else:
                raise RuntimeError(f"Keyword SRCTYPE is {srctype}.  It should "
                    "be 'POINT' or 'EXTENDED'. Can't decide between `flux` and "
                    "`surf_bright` columns.")

            if np.min(uncertainty.array) <= 0.:
                warnings.warn("Standard Deviation has values of 0 or less",
                    AstropyUserWarning)

            meta = dict(slitname=hdu.header.get('SLTNAME', ''))

            spec = Spectrum1D(flux=flux, spectral_axis=wavelength,
                uncertainty=uncertainty, meta=meta)
            spectra.append(spec)

    return SpectrumList(spectra)
