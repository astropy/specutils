import warnings

import astropy.units as u
from astropy.units import Quantity
from astropy.table import Table
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.utils.exceptions import AstropyUserWarning
import numpy as np
import asdf
from gwcs.wcstools import grid_from_bounding_box

from ...spectra import Spectrum1D, SpectrumList
from ..registers import data_loader


__all__ = ["jwst_x1d_single_loader", "jwst_x1d_multi_loader"]


def identify_jwst_x1d_fits(origin, *args, **kwargs):
    """
    Check whether the given file is a JWST x1d spectral data product.
    """
    is_jwst = _identify_jwst_fits(args[0])
    with fits.open(args[0], memmap=False) as hdulist:
        if (is_jwst and 'EXTRACT1D' in hdulist and ('EXTRACT1D', 2) not in hdulist
            and "SCI" not in hdulist):
            return True
        else:
            return False


def identify_jwst_x1d_multi_fits(origin, *args, **kwargs):
    """
    Check whether the given file is a JWST x1d spectral data product.
    """
    is_jwst = _identify_jwst_fits(args[0])
    with fits.open(args[0], memmap=False) as hdulist:
        if is_jwst and ('EXTRACT1D', 2) in hdulist and "SCI" not in hdulist:
            return True
        else:
            return False


def identify_jwst_s2d_fits(origin, *args, **kwargs):
    """
    Check whether the given file is a JWST s2d spectral data product.
    """
    is_jwst = _identify_jwst_fits(args[0])
    with fits.open(args[0], memmap=False) as hdulist:
        if (is_jwst and "SCI" in hdulist and ("SCI", 2) not in hdulist
            and "EXTRACT1D" not in hdulist):
            return True
        else:
            return False


def identify_jwst_s2d_multi_fits(origin, *args, **kwargs):
    """
    Check whether the given file is a JWST s2d spectral data product.
    """
    is_jwst = _identify_jwst_fits(args[0])
    with fits.open(args[0], memmap=False) as hdulist:
        if is_jwst and ("SCI", 2) in hdulist and "EXTRACT1D" not in hdulist:
            return True
        else:
            return False


def _identify_jwst_fits(filename):
    """
    Check whether the given file is a JWST spectral data product.
    """
    try:
        with fits.open(filename, memmap=False) as hdulist:
            # This is a near-guarantee that we have a JWST data product
            if 'ASDF' not in hdulist:
                return False
            if not hdulist[0].header["TELESCOP"] == ("JWST"):
                return False
        return True
    # This probably means we didn't have a FITS file
    except Exception:
        return False


@data_loader("JWST x1d", identifier=identify_jwst_x1d_fits, dtype=Spectrum1D,
            extensions=['fits'])
def jwst_x1d_single_loader(filename, **kwargs):
    """
    Loader for JWST x1d 1-D spectral data in FITS format

    Parameters
    ----------
    filename: str
        The path to the FITS file

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file.
    """
    spectrum_list = _jwst_x1d_loader(filename, **kwargs)
    if len(spectrum_list) == 1:
        return spectrum_list[0]
    else:
        raise RuntimeError(f"Input data has {len(spectrum_list)} spectra. "
            "Use SpectrumList.read() instead.")


@data_loader("JWST x1d multi", identifier=identify_jwst_x1d_multi_fits,
            dtype=SpectrumList, extensions=['fits'])
def jwst_x1d_multi_loader(filename, **kwargs):
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
    return _jwst_x1d_loader(filename, **kwargs)


def _jwst_x1d_loader(filename, **kwargs):
    """Implementation of loader for JWST x1d 1-D spectral data in FITS format

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

    with fits.open(filename, memmap=False) as hdulist:

        primary_header = hdulist[0].header

        for hdu in hdulist:
            # Read only the BinaryTableHDUs named EXTRACT1D
            if hdu.name != 'EXTRACT1D':
                continue

            data = Table.read(hdu)

            wavelength = Quantity(data["WAVELENGTH"])

            # Determine if FLUX or SURF_BRIGHT column should be returned
            # based on whether it is point or extended source
            try:
                srctype = hdu.header.get("srctype")
            except KeyError:
                # SRCTYPE is in primary header because there is only one spectrum
                srctype = hdulist.header.get("srctype")

            if srctype == "POINT":
                flux = Quantity(data["FLUX"])
                uncertainty = StdDevUncertainty(data["ERROR"])

            elif srctype == "EXTENDED":
                flux = Quantity(data["SURF_BRIGHT"])
                uncertainty = StdDevUncertainty(hdu.data["SB_ERROR"])

            else:
                raise RuntimeError(f"Keyword SRCTYPE is {srctype}.  It should "
                    "be 'POINT' or 'EXTENDED'. Can't decide between `flux` and "
                    "`surf_bright` columns.")

            # Merge primary and slit headers and dump into meta
            slit_header = hdu.header
            header = primary_header.copy()
            header.extend(slit_header, strip=True, update=True)
            meta = {k: v for k,v in header.items()}

            spec = Spectrum1D(flux=flux, spectral_axis=wavelength,
                uncertainty=uncertainty, meta=meta)
            spectra.append(spec)

    return SpectrumList(spectra)


@data_loader("JWST s2d", identifier=identify_jwst_s2d_fits, dtype=Spectrum1D,
            extensions=['fits'])
def jwst_s2d_single_loader(filename, **kwargs):
    """
    Loader for JWST s2d 2D rectified spectral data in FITS format

    Parameters
    ----------
    filename: str
        The path to the FITS file

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file.
    """
    spectrum_list = _jwst_s2d_loader(filename, **kwargs)
    if len(spectrum_list) == 1:
        return spectrum_list[0]
    elif len(spectrum_list) > 1:
        raise RuntimeError(f"Input data has {len(spectrum_list)} spectra. "
            "Use SpectrumList.read() instead.")
    else:
        raise RuntimeError(f"Input data has {len(spectrum_list)} spectra.")


@data_loader("JWST s2d multi", identifier=identify_jwst_s2d_multi_fits, dtype=SpectrumList,
            extensions=['fits'])
def jwst_s2d_multi_loader(filename, **kwargs):
    """
    Loader for JWST s2d 2D rectified spectral data in FITS format

    Parameters
    ----------
    filename: str
        The path to the FITS file

    Returns
    -------
    SpectrumList
        The spectra contained in the file.
    """
    return _jwst_s2d_loader(filename, **kwargs)


def _jwst_s2d_loader(filename, **kwargs):
    """
    Loader for JWST s2d 2D rectified spectral data in FITS format

    Parameters
    ----------
    filename: str
        The path to the FITS file

    Returns
    -------
    SpectrumList
        The spectra contained in the file.
    """
    spectra = []

    # Get a list of GWCS objects from the slits
    with asdf.open(filename) as af:
        # Slits can be listed under "slits", "products" or "exposures"
        if "products" in af.tree:
            slits = "products"
        elif "exposures" in af.tree:
            slits = "exposures"
        elif "slits" in af.tree:
            slits = "slits"
        else:
            raise RuntimeError(f"Cannot load gwcs object in {filename}")

        # Create list of the GWCS objects, one for each slit
        wcslist = [slit["meta"]["wcs"] for slit in af.tree[slits]]

    with fits.open(filename, memmap=False) as hdulist:

        primary_header = hdulist["PRIMARY"].header

        hdulist_sci = [hdu for hdu in hdulist if hdu.name == "SCI"]

        for hdu, wcs in zip(hdulist_sci, wcslist):
            # Get flux
            try:
                flux_unit = u.Unit(hdu.header["BUNIT"])
            except (ValueError, KeyError):
                flux_unit = u.Unit("MJy")
            flux = Quantity(hdu.data, unit=flux_unit)

            # Get the wavelength array from the GWCS object which returns a
            # tuple of (RA, Dec, lambda)
            grid = grid_from_bounding_box(wcs.bounding_box)
            _, _, lam = wcs(*grid)
            _, _, lam_unit = wcs.output_frame.unit

            # The dispersion axis is 1-indexed in the FITS file, but we
            # need it zero-indexed
            dispaxis = hdu.header["DISPAXIS"] - 1

            # Make sure the dispersion axis is the same for every spatial axis
            if not (lam == lam[dispaxis]).all():
                raise RuntimeError("This 2D or 3D spectrum is not rectified "
                    "and cannot be loaded into a Spectrum1D object")
            wavelength = Quantity(lam[dispaxis], unit=lam_unit)

            # Merge primary and slit headers and dump into meta
            slit_header = hdu.header
            header = primary_header.copy()
            header.extend(slit_header, strip=True, update=True)
            meta = {k: v for k,v in header.items()}

            spec = Spectrum1D(flux=flux, spectral_axis=wavelength, meta=meta)
            spectra.append(spec)

    return SpectrumList(spectra)
