"""Register reader functions for various spectral formats."""
from collections import OrderedDict
from functools import wraps
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.nddata import StdDevUncertainty, InverseVariance
from specutils import SpectralAxis, Spectrum1D, SpectrumList
from specutils.io.registers import get_loaders_by_extension, io_registry
"""
From the specutils documentation:
- For spectra that have different shapes, use SpectrumList.
- For spectra that have the same shape but different spectral axes, see SpectrumCollection. 
- For a spectrum or spectra that all share the same spectral axis, use Spectrum1D. 
"""


## CLEANING REGISTRY
def clear_fits_registry_for_spectrum_objects(extension, dtypes):
    for data_format in get_loaders_by_extension(extension):
        for dtype in dtypes:
            try:
                io_registry.unregister_identifier(data_format, dtype)
            except:
                continue
    return None


# The registry is full of junk, and many of them don't even work out of the box.
clear_fits_registry_for_spectrum_objects(extension="fits",
                                         dtypes=(Spectrum1D, SpectrumList))


## DATA LOADERS
def data_loader(label,
                identifier,
                dtype,
                extensions=None,
                priority=0,
                force=False):

    def identifier_wrapper(ident):

        def wrapper(*args, **kwargs):
            try:
                return ident(*args, **kwargs)
            except Exception as e:
                return False

        return wrapper

    def decorator(func):
        io_registry.register_reader(label,
                                    dtype,
                                    func,
                                    priority=priority,
                                    force=force)
        io_registry.register_identifier(label,
                                        dtype,
                                        identifier_wrapper(identifier),
                                        force=force)
        func.extensions = extensions

        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)

        return wrapper

    return decorator


## HELPERS

# TODO: rewrite into actual function
is_filetype = lambda filetype: lambda f, o, *args, **kwargs: o.split("/")[
    -1].startswith(filetype)


def _wcs_log_linear(naxis, cdelt, crval):
    """
    Convert WCS from log to linear.

    Note: WCS used is SDSS format, not standard FITS format.
    """
    return 10**(np.arange(naxis) * cdelt + crval)


def _fetch_metadata(image: fits.hdu.hdulist.HDUList):
    """
    Fetch the relevant metadata.

    Global function used in all SDSS-V metadata handling.

    Parameters
    ----------
    image

    Returns
    -------
    meta

    """
    # Orderly with an OrderedDict
    common_meta = OrderedDict([])

    # Loop over all keys to get relevant metadata
    for key in image[0].header.keys():
        # exclude these keys
        if key.startswith(("TTYPE", "TFORM", "TDIM")) or key in (
                "",
                "COMMENT",
                "CHECKSUM",
                "DATASUM",
                "NAXIS",
                "NAXIS1",
                "NAXIS2",
                "XTENSION",
                "BITPIX",
                "PCOUNT",
                "GCOUNT",
                "TFIELDS",
        ):
            continue

        common_meta[key.lower()] = image[0].header[key]  # add key to dict

    return common_meta


def _fetch_flux_unit(hdu: fits.hdu.image.ImageHDU):
    """
    Fetch the flux unit by accessing the BUNIT key of the given HDU.

    Parameters
    ----------
    hdu : flux or e_flux hdu

    Returns
    -------
    flux_unit : astropy Units object of flux unit

    """
    flux_unit = (hdu.header["BUNIT"].replace("Ang", "Angstrom").replace(
        "(", "").replace(")", "").split(" ")[-2:])

    # TODO: make it so the stupid warning is surpressed
    """

    """
    return u.Unit("".join(flux_unit))


## APOGEE files
# TODO: add LSF and wavelength coefficients to metadata?
# @data_loader(
#    "apStar",
#    identifier=is_filetype(("apStar", "asStar")),
#    dtype=Spectrum1D,
#    priority=10,
#    extensions=["fits"],
# )
def load_sdss_apStar_1D(path, data_slice=None, **kwargs):
    # intialize slicer
    slicer = slice(*data_slice) if data_slice is not None else slice(None)

    with fits.open(path) as image:
        # Obtain spectral axis from header (HDU[1])
        wavelength = _wcs_log_linear(
            image[1].header["NAXIS1"],
            image[1].header["CDELT1"],
            image[1].header["CRVAL1"],
        )
        spectral_axis = u.Quantity(wavelength, unit=u.Angstrom)

        # Obtain flux and e_flux from header (HDU[1] & HDU[2])
        flux_unit = _fetch_flux_unit(image[1])
        flux = u.Quantity(image[1].data[slicer], unit=flux_unit)
        e_flux = StdDevUncertainty(image[2].data[slicer])

        snr = [image[0].header["SNR"]]
        n_visits = image[0].header["NVISITS"]
        if n_visits > 1:
            snr.append(
                snr[0])  # duplicate S/N value for second stacking method
            snr.extend([
                image[0].header[f"SNRVIS{i}"] for i in range(1, 1 + n_visits)
            ])

        meta = _fetch_metadata(image)

        meta["SNR"] = np.array(snr)[slicer]
        meta["BITMASK"] = image[3].data[slicer]

    return Spectrum1D(spectral_axis=spectral_axis,
                      flux=flux,
                      uncertainty=e_flux,
                      meta=meta)

    # @data_loader(
    #    "apStar",
    #    identifier=is_filetype(("apStar", "asStar")),
    #    dtype=SpectrumList,
    #    priority=10,
    #    extensions=["fits"],
    # )


def load_sdss_apStar_list(path, **kwargs):
    return SpectrumList([load_sdss_apStar(path, **kwargs)])


# @data_loader(
#    "apVisit",
#    identifier=is_filetype("apVisit"),
#    dtype=Spectrum1D,
#    priority=10,
#    extensions=["fits"],
# )
def load_sdss_apVisit_1D(path, **kwargs):
    flux_unit = u.Unit("1e-17 erg / (Angstrom cm2 s)")

    with fits.open(path) as image:
        flux_unit = _fetch_flux_unit(image[1])
        # create sorted data array
        spectral_axis = u.Quantity(image[4].data.flatten(), unit=u.Angstrom)
        print(spectral_axis)

        # Handle chips
        flux = u.Quantity(image[1].data.flatten(), unit=flux_unit)
        e_flux = StdDevUncertainty(image[2].data.flatten())

        # Get metadata
        meta = _fetch_metadata(image)

        meta["bitmask"] = image[3].data.flatten()
        # TODO: Include things like sky flux, sky error, telluric flux, telluric error?
        #       wavelength coefficients? lsf coefficients?

    return Spectrum1D(spectral_axis=spectral_axis,
                      flux=flux,
                      uncertainty=e_flux,
                      meta=meta)

    # @data_loader(
    #    "apVisit",
    #    identifier=is_filetype("apVisit"),
    #    dtype=SpectrumList,
    #    priority=10,
    #    extensions=["fits"],
    # )


def load_sdss_apVisit_multi(path, **kwargs):
    spectra = SpectrumList()
    with fits.open(path) as image:
        # Get metadata
        flux_unit = _fetch_flux_unit(image[1])
        common_meta = _fetch_metadata(image)

        for chip in range(image[1].data.shape[0]):
            spectral_axis = u.Quantity(image[4].data[chip], unit=u.Angstrom)
            flux = u.Quantity(image[1].data[chip], unit=flux_unit)
            e_flux = StdDevUncertainty(image[2].data[chip])

            meta = common_meta.copy()
            meta["BITMASK"] = image[3].data[chip]

            # TODO: Include things like sky flux, sky error, telluric flux, telluric error?
            #       wavelength coefficients? lsf coefficients?

            spectra.append(
                Spectrum1D(
                    spectral_axis=spectral_axis,
                    flux=flux,
                    uncertainty=e_flux,
                    meta=meta,
                ))

    return spectra


## BOSS files (specLite, specFull)


# @data_loader(
#    "specFull",
#    # Note the path definition here is not the same as other SDSS-V data models.
#    identifier=is_filetype("spec"),
#    dtype=Spectrum1D,
#    priority=10,
#    extensions=["fits"],
# )
def load_sdss_specFull_1D(path, **kwargs):
    with fits.open(path) as image:
        # Fetch flux, e_flux, spectral axis, and units
        # TODO: check with data team about loc of BUNIT in specFull/Lite files
        flux_unit = u.Unit("1e-17 erg / (Angstrom cm2 s)")
        spectral_axis = u.Quantity(10**image[1].data["LOGLAM"],
                                   unit=u.Angstrom)

        flux = u.Quantity(image[1].data["FLUX"], unit=flux_unit)
        # no e_flux, so we use inverse of variance
        ivar = InverseVariance(image[1].data["IVAR"])

        # Fetch metadata
        # Note: specFull file does not include S/N value, but this gets calculated
        #       for mwmVisit/mwmStar files when they are created
        meta = _fetch_metadata(image)

    return Spectrum1D(spectral_axis=spectral_axis,
                      flux=flux,
                      uncertainty=ivar,
                      meta=meta)

    # @data_loader(
    #    "specFull",
    #    # Note the path definition here is not the same as other SDSS-V data models.
    #    identifier=is_filetype("spec"),
    #    dtype=SpectrumList,
    #    priority=10,
    #    extensions=["fits"],
    # )


def load_sdss_specFull_list(path, **kwargs):
    with fits.open(path) as image:
        spectra = SpectrumList()
        for i in range(len(image) - 1):
            # skip non-exposures
            # TODO: I have no idea about these formats, need to check what's up with them
            if image[i + 1].name in ("SPALL", "ZALL", "ZLINE"):
                continue

            # Fetch flux, e_flux, spectral axis, and units
            # TODO: check with data team about loc of BUNIT in specFull/Lite files

            flux_unit = u.Unit("1e-17 erg / (Angstrom cm2 s)")
            spectral_axis = u.Quantity(10**image[i + 1].data["LOGLAM"],
                                       unit=u.Angstrom)

            flux = u.Quantity(image[i + 1].data["FLUX"], unit=flux_unit)
            # no e_flux, so we use inverse of variance
            ivar = InverseVariance(image[i + 1].data["IVAR"])

            # Fetch metadata
            # Note: specFull file does not include S/N value, but this gets calculated
            #       for mwmVisit/mwmStar files when they are created
            meta = _fetch_metadata(image)
            meta["name"] = image[i + 1].name

            spectra.append(
                Spectrum1D(spectral_axis=spectral_axis,
                           flux=flux,
                           uncertainty=ivar,
                           meta=meta))

    return spectra
