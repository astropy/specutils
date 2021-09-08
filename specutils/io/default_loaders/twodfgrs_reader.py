from astropy.nddata import VarianceUncertainty
from astropy.units import Unit
from astropy.wcs import WCS

from ...spectra import Spectrum1D, SpectrumList
from ..registers import data_loader
from ..parsing_utils import read_fileobj_or_hdulist


def identify_2dfgrs(origin, *args, **kwargs):
    """
    Identify if the current file is a 2dFGRS file
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return hdulist[0].header.get("IMAGE", "").strip() in ("SKYCHART", "SPECTRUM")


def load_spectrum_from_extension(hdu, primary_header):
    """
    2dFGRS files contain multiple spectra, this is a helper function to load the
    spectrum for a single HDU.
    """
    header = hdu.header
    # Due to the missing information, the WCS needs to be read in manually
    wcs = WCS(naxis=1)
    wcs.wcs.cdelt[0] = header["CDELT1"]
    wcs.wcs.crval[0] = header["CRVAL1"]
    wcs.wcs.crpix[0] = header["CRPIX1"]
    wcs.wcs.cunit[0] = Unit("Angstrom")

    meta = {
        "header": header,
        "primary_header": primary_header,
    }

    spectrum = hdu.data[0] * Unit("count/s")
    variance = VarianceUncertainty(hdu.data[1])
    sky = hdu.data[2] * Unit("count/s")

    meta["sky"] = Spectrum1D(flux=sky, wcs=wcs)

    return Spectrum1D(flux=spectrum, wcs=wcs, meta=meta, uncertainty=variance)


@data_loader(
    "2dFGRS", identifier=identify_2dfgrs, dtype=SpectrumList,
    extensions=["fit", "fits"], priority=10,
)
def twodfgrs_fits_loader(file_obj, **kwargs):
    """
    Load a file from the 2dF Galaxy Redshift Survey.

    The 2dF Galaxy Redshift Survey (2dFGRS) was a major spectroscopic survey
    taking full advantage of the unique capabilities of the 2dF facility built
    by the Anglo-Australian Observatory. The 2dFGRS is integrated with the 2dF
    QSO survey. Further details can be seen at http://www.2dfgrs.net/ or at
    https://docs.datacentral.org.au/.

    Parameters
    ----------
    file_obj: str or HDUList
        The path to the FITS file, or the loaded file as an HDUList
    Returns
    -------
    data: SpectrumList
        The 2dFGRS spectra represented by the data in this file.
    """

    spectra = []

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        primary_header = hdulist[0].header

        for hdu in hdulist:
            if hdu.header.get("IMAGE", "").strip() == "SKYCHART":
                continue
            spectra.append(load_spectrum_from_extension(hdu, primary_header))

    return SpectrumList(spectra)
