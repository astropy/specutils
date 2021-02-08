import astropy.units as u
from astropy.units import Quantity
from astropy.table import Table
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
import numpy as np
import asdf
from gwcs.wcstools import grid_from_bounding_box

from ...spectra import Spectrum1D, SpectrumList
from ..registers import data_loader
from ..parsing_utils import read_fileobj_or_hdulist


__all__ = ["jwst_x1d_single_loader", "jwst_x1d_multi_loader"]


def identify_jwst_x1d_fits(origin, *args, **kwargs):
    """
    Check whether the given file is a JWST x1d spectral data product.
    """
    is_jwst = _identify_jwst_fits(*args)
    with read_fileobj_or_hdulist(*args, memmap=False, **kwargs) as hdulist:
        return (is_jwst and 'EXTRACT1D' in hdulist and ('EXTRACT1D', 2) not in hdulist)


def identify_jwst_x1d_multi_fits(origin, *args, **kwargs):
    """
    Check whether the given file is a JWST x1d spectral data product with many slits.
    """
    is_jwst = _identify_jwst_fits(*args)
    with read_fileobj_or_hdulist(*args, memmap=False, **kwargs) as hdulist:
        return is_jwst and ('EXTRACT1D', 2) in hdulist


def identify_jwst_s2d_fits(origin, *args, **kwargs):
    """
    Check whether the given file is a JWST s2d spectral data product.
    """
    is_jwst = _identify_jwst_fits(*args)
    with read_fileobj_or_hdulist(*args, memmap=False, **kwargs) as hdulist:
        return (is_jwst and "SCI" in hdulist and ("SCI", 2) not in hdulist
                and "EXTRACT1D" not in hdulist and len(hdulist["SCI"].data.shape) == 2)


def identify_jwst_s2d_multi_fits(origin, *args, **kwargs):
    """
    Check whether the given file is a JWST s2d spectral data product with many slits.
    """
    is_jwst = _identify_jwst_fits(*args)
    with read_fileobj_or_hdulist(*args, memmap=False, **kwargs) as hdulist:
        return (is_jwst and ("SCI", 2) in hdulist and "EXTRACT1D" not in hdulist
                and len(hdulist["SCI"].data.shape) == 2)


def identify_jwst_s3d_fits(origin, *args, **kwargs):
    """
    Check whether the given file is a JWST s3d spectral data product.
    """
    is_jwst = _identify_jwst_fits(*args)
    with read_fileobj_or_hdulist(*args, memmap=False, **kwargs) as hdulist:
        return (is_jwst and "SCI" in hdulist and "EXTRACT1D" not in hdulist
                and len(hdulist["SCI"].data.shape) == 3)


def _identify_jwst_fits(*args):
    """
    Check whether the given file is a JWST data product.
    """
    try:
        with read_fileobj_or_hdulist(*args, memmap=False) as hdulist:
            return "ASDF" in hdulist and hdulist[0].header.get("TELESCOP") == "JWST"
    # This probably means we didn't have a FITS file
    except Exception:
        return False


@data_loader("JWST x1d", identifier=identify_jwst_x1d_fits, dtype=Spectrum1D,
             extensions=['fits'])
def jwst_x1d_single_loader(file_obj, **kwargs):
    """
    Loader for JWST x1d 1-D spectral data in FITS format

    Parameters
    ----------
    filename : str
        The path to the FITS file

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file.
    """
    spectrum_list = _jwst_x1d_loader(file_obj, **kwargs)
    if len(spectrum_list) == 1:
        return spectrum_list[0]
    else:
        raise RuntimeError(f"Input data has {len(spectrum_list)} spectra. "
                           "Use SpectrumList.read() instead.")


@data_loader("JWST x1d multi", identifier=identify_jwst_x1d_multi_fits,
             dtype=SpectrumList, extensions=['fits'])
def jwst_x1d_multi_loader(file_obj, **kwargs):
    """
    Loader for JWST x1d 1-D spectral data in FITS format

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
          FITS file name, object (provided from name by Astropy I/O Registry),
          or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    SpectrumList
        A list of the spectra that are contained in the file.
    """
    return _jwst_x1d_loader(file_obj, **kwargs)


def _jwst_x1d_loader(file_obj, **kwargs):
    """Implementation of loader for JWST x1d 1-D spectral data in FITS format

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
          FITS file name, object (provided from name by Astropy I/O Registry),
          or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    SpectrumList
        A list of the spectra that are contained in the file.
    """

    spectra = []

    with read_fileobj_or_hdulist(file_obj, memmap=False, **kwargs) as hdulist:

        primary_header = hdulist["PRIMARY"].header

        for hdu in hdulist:
            # Read only the BinaryTableHDUs named EXTRACT1D and SCI
            if hdu.name != 'EXTRACT1D':
                continue

            data = Table.read(hdu)

            wavelength = Quantity(data["WAVELENGTH"])

            # Determine if FLUX or SURF_BRIGHT column should be returned
            # based on whether it is point or extended source.
            #
            # SRCTYPE used to be in primary header, but was moved around. As
            # per most recent pipeline definition, it should be in the
            # EXTRACT1D extension.
            srctype = None
            if "srctype" in hdu.header:
                srctype = hdu.header.get("srctype")

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
            meta = {'header': header}

            spec = Spectrum1D(flux=flux, spectral_axis=wavelength,
                              uncertainty=uncertainty, meta=meta)
            spectra.append(spec)

    return SpectrumList(spectra)


@data_loader("JWST s2d", identifier=identify_jwst_s2d_fits, dtype=Spectrum1D,
             extensions=['fits'])
def jwst_s2d_single_loader(filename, **kwargs):
    """
    Loader for JWST s2d 2D rectified spectral data in FITS format.

    Parameters
    ----------
    filename : str
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
    Loader for JWST s2d 2D rectified spectral data in FITS format.

    Parameters
    ----------
    filename : str
        The path to the FITS file

    Returns
    -------
    SpectrumList
        The spectra contained in the file.
    """
    return _jwst_s2d_loader(filename, **kwargs)


def _jwst_s2d_loader(filename, **kwargs):
    """
    Loader for JWST s2d 2D rectified spectral data in FITS format.

    Parameters
    ----------
    filename : str
        The path to the FITS file

    Returns
    -------
    SpectrumList
        The spectra contained in the file.
    """
    spectra = []
    slits = None

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
            # No slits.  Case for s3d cubes
            wcslist = [af.tree["meta"]["wcs"]]

        # Create list of the GWCS objects, one for each slit
        if slits is not None:
            wcslist = [slit["meta"]["wcs"] for slit in af.tree[slits]]

    with fits.open(filename, memmap=False) as hdulist:

        primary_header = hdulist["PRIMARY"].header

        hdulist_sci = [hdu for hdu in hdulist if hdu.name == "SCI"]

        for hdu, wcs in zip(hdulist_sci, wcslist):
            # Get flux
            try:
                flux_unit = u.Unit(hdu.header["BUNIT"])
            except (ValueError, KeyError):
                flux_unit = None

            # The dispersion axis is 1 or 2.  1=x, 2=y.
            dispaxis = hdu.header.get("DISPAXIS")
            if dispaxis is None:
                dispaxis = hdulist["PRIMARY"].header.get("DISPAXIS")

            # Get the wavelength array from the GWCS object which returns a
            # tuple of (RA, Dec, lambda)
            grid = grid_from_bounding_box(wcs.bounding_box)
            _, _, lam = wcs(*grid)
            _, _, lam_unit = wcs.output_frame.unit

            # Make sure the dispersion axis is the same for every spatial axis
            # for s2d data
            if dispaxis == 1:
                flux_array = hdu.data
                wavelength_array = lam[0]
                # Make sure all rows are the same
                if not (lam == wavelength_array).all():
                    raise RuntimeError("This 2D or 3D spectrum is not rectified "
                                       "and cannot be loaded into a Spectrum1D object.")
            elif dispaxis == 2:
                flux_array = hdu.data.T
                wavelength_array = lam[:, 0]
                # Make sure all columns are the same
                if not (lam.T == lam[None, :, 0]).all():
                    raise RuntimeError("This 2D or 3D spectrum is not rectified "
                                       "and cannot be loaded into a Spectrum1D object.")
            else:
                raise RuntimeError("This 2D spectrum has an unknown dispaxis "
                                   "and cannot be loaded into a Spectrum1D object.")

            flux = Quantity(flux_array, unit=flux_unit)
            wavelength = Quantity(wavelength_array, unit=lam_unit)

            # Merge primary and slit headers and dump into meta
            slit_header = hdu.header
            header = primary_header.copy()
            header.extend(slit_header, strip=True, update=True)
            meta = {'header': header}

            spec = Spectrum1D(flux=flux, spectral_axis=wavelength, meta=meta)
            spectra.append(spec)

    return SpectrumList(spectra)


@data_loader("JWST s3d", identifier=identify_jwst_s3d_fits, dtype=Spectrum1D,
             extensions=['fits'])
def jwst_s3d_single_loader(filename, **kwargs):
    """
    Loader for JWST s3d 3D rectified spectral data in FITS format.

    Parameters
    ----------
    filename : str
        The path to the FITS file

    Returns
    -------
    Spectrum1D
        The spectrum contained in the file.
    """
    spectrum_list = _jwst_s3d_loader(filename, **kwargs)
    if len(spectrum_list) == 1:
        return spectrum_list[0]
    elif len(spectrum_list) > 1:
        raise RuntimeError(f"Input data has {len(spectrum_list)} spectra. "
                           "Use SpectrumList.read() instead.")
    else:
        raise RuntimeError(f"Input data has {len(spectrum_list)} spectra.")


def _jwst_s3d_loader(filename, **kwargs):
    """
    Loader for JWST s3d 3D rectified spectral data in FITS format.

    Parameters
    ----------
    filename : str
        The path to the FITS file

    Returns
    -------
    SpectrumList
        The spectra contained in the file.
    """
    spectra = []

    # Get a list of GWCS objects from the slits
    with asdf.open(filename) as af:
        wcslist = [af.tree["meta"]["wcs"]]

    with fits.open(filename, memmap=False) as hdulist:

        primary_header = hdulist["PRIMARY"].header

        hdulist_sci = [hdu for hdu in hdulist if hdu.name == "SCI"]

        for hdu, wcs in zip(hdulist_sci, wcslist):
            # Get flux
            try:
                flux_unit = u.Unit(hdu.header["BUNIT"])
            except (ValueError, KeyError):
                flux_unit = None

            # The spectral axis is first.  We need it last
            flux_array = hdu.data.T
            flux = Quantity(flux_array, unit=flux_unit)

            # Get the wavelength array from the GWCS object which returns a
            # tuple of (RA, Dec, lambda)
            grid = grid_from_bounding_box(wcs.bounding_box)
            _, _, lam = wcs(*grid)
            _, _, lam_unit = wcs.output_frame.unit

            wavelength_array = lam[:, 0, 0]
            wavelength = Quantity(wavelength_array, unit=lam_unit)

            # Merge primary and slit headers and dump into meta
            slit_header = hdu.header
            header = primary_header.copy()
            header.extend(slit_header, strip=True, update=True)
            meta = {'header': header}

            spec = Spectrum1D(flux=flux, spectral_axis=wavelength, meta=meta)
            spectra.append(spec)

    return SpectrumList(spectra)
