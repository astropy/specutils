from astropy.units import Unit
from astropy.wcs import WCS

from ...spectra import Spectrum1D, SpectrumList
from ..registers import data_loader
from ..parsing_utils import read_fileobj_or_hdulist


def identify_2slaq_lrg(origin, *args, **kwargs):
    """
    Identify if the current file is a 2SLAQ-LRG file
    """
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        if hdulist[0].header["MSTITLE"].startswith("2dF-SDSS LRG/QSO survey"):
            # apparently the QSO part doesn't have MSTITLE, so we should be safe
            # with just the above condition, but just in case, we know they have
            # a different structure (LRG has one ext, QSO has more; LRG is 3d,
            # QSO is 1d
            if len(hdulist) == 1 and hdulist[0].data.ndim == 3:
                return True


@data_loader(
    "2SLAQ-LRG", identifier=identify_2slaq_lrg, dtype=SpectrumList,
    extensions=["fit", "fits"], priority=10,
)
def twoslaq_lrg_fits_loader(file_obj, **kwargs):
    """
    Load a file from the LRG subset of the 2dF-SDSS LRG/QSO survey (2SLAQ-LRG)
    file.

    2SLAQ was one of a number of surveys that used the 2dF instrument on the
    Anglo-Australian Telescope (AAT) at Siding Spring Observatory (SSO) near
    Coonabarabran, Australia, and followed up the 2QZ survey. Further details
    can be seen at http://www.physics.usyd.edu.au/2slaq/ or at
    https://docs.datacentral.org.au/.

    The LRG and QSO data appear to be in different formats, this only loads the
    LRG version. As there is a science and sky spectrum, the `SpectrumList`
    loader provides both, whereas the `Spectrum1D` loader only provides the
    science.

    Parameters
    ----------
    file_name: str
        The path to the FITS file
    Returns
    -------
    data: SpectrumList
        The 2SLAQ-LRG spectrum that is represented by the data in this file.
    """

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        header = hdulist[0].header
        spectrum = hdulist[0].data[0, 0] * Unit("count/s")
        sky = hdulist[0].data[1, 0] * Unit("count/s")

    # Due to the odd structure of the file, the WCS needs to be read in
    # manually
    wcs = WCS(naxis=1)
    wcs.wcs.cdelt[0] = header["CD1_1"]
    wcs.wcs.crval[0] = header["CRVAL1"]
    wcs.wcs.crpix[0] = header["CRPIX1"]
    wcs.wcs.cunit[0] = Unit("Angstrom")

    meta = {"header": header}

    return SpectrumList([
        Spectrum1D(flux=spectrum, wcs=wcs, meta=meta),
        Spectrum1D(flux=sky, wcs=wcs, meta=meta),
    ])

# Commented out until discussion about whether to provide science-only or not
# @data_loader("2SLAQ-LRG", identifier=identify_2slaq_lrg,dtype=Spectrum1D,
#              extensions=["fit", "fits"])
# def twoslaq_lrg_fits_loader_only_science(filename, **kwargs):
#     """
#     Load a file from the LRG subset of the 2dF-SDSS LRG/QSO survey (2SLAQ-LRG)
#     file.
#
#     2SLAQ was one of a number of surveys that used the 2dF instrument on the
#     Anglo-Australian Telescope (AAT) at Siding Spring Observatory (SSO) near
#     Coonabarabran, Australia, and followed up the 2QZ survey. Further details
#     can be seen at http://www.physics.usyd.edu.au/2slaq/ or at
#     https://docs.datacentral.org.au/.
#
#     The LRG and QSO data appear to be in different formats, this only loads the
#     LRG version. As there is a science and sky spectrum, the `SpectrumList`
#     loader provides both, whereas the `Spectrum1D` loader only provides the
#     science.
#
#     Parameters
#     ----------
#     file_name: str
#         The path to the FITS file
#     Returns
#     -------
#     data: Spectrum1D
#         The 2SLAQ-LRG spectrum that is represented by the data in this file.
#     """
#     return SpectrumList.read(filename, format="2SLAQ-LRG")[0]
