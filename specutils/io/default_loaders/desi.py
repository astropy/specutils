"""
Loader for DESI spectrum files: spectra_ and coadd_.

* spectra_ files contain all observations of an object in a given region.
* coadd_ files contain one, co-added spectrum of an object in a given region.
  The coaddition is performed across observations, but *not* across different
  cameras.

.. _spectra: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/spectra-SURVEY-PROGRAM-PIXNUM.html
.. _coadd: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/coadd-SURVEY-PROGRAM-PIXNUM.html
"""
import re

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.units import Unit
from astropy.nddata import StdDevUncertainty, InverseVariance

import numpy as np

from ...spectra import Spectrum1D, SpectrumList
from ..registers import data_loader
from ..parsing_utils import read_fileobj_or_hdulist

__all__ = ['spectra_identify', 'coadd_identify'
           'spectra_loader', 'coadd_loader']

#
# These are reserved for future use.
#
_spectra_pattern = re.compile(r'spectra-(cmx|main|special|sv1|sv2|sv3)-(backup|bright|dark|other)-[0-9]+\.fits')
_coadd_pattern = re.compile(r'coadd-(cmx|main|special|sv1|sv2|sv3)-(backup|bright|dark|other)-[0-9]+\.fits')


def spectra_identify(origin, *args, **kwargs):
    """
    Check whether given input is FITS and has DESI stuff.
    This is used for Astropy I/O Registry.
    """
    # Test if fits has extension of type BinTable and check for spec-specific keys
    # with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
    #     return (hdulist[0].header.get('TELESCOP') == 'SDSS 2.5-M' and
    #             hdulist[0].header.get('FIBERID', 0) > 0 and
    #             len(hdulist) > 1 and
    #             (isinstance(hdulist[1], fits.BinTableHDU) and
    #              hdulist[1].header.get('TTYPE3').lower() == 'ivar'))
    return True


def coadd_identify(origin, *args, **kwargs):
    """
    Check whether given input is FITS and has DESI stuff.
    This is used for Astropy I/O Registry.
    """
    # Test if fits has extension of type BinTable and check for spec-specific keys
    # with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
    #     return (hdulist[0].header.get('TELESCOP') == 'SDSS 2.5-M' and
    #             hdulist[0].header.get('FIBERID', 0) > 0 and
    #             len(hdulist) > 1 and
    #             (isinstance(hdulist[1], fits.BinTableHDU) and
    #              hdulist[1].header.get('TTYPE3').lower() == 'ivar'))
    return True


@data_loader(
    label="DESI spectra", identifier=spectra_identify, extensions=['fits'],
    priority=10,
)
def spectra_loader(file_obj, **kwargs):
    """
    Loader for DESI spectra_ files.

    .. _spectra: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/spectra-SURVEY-PROGRAM-PIXNUM.html

    Any keyword arguments are passed to
    :func:`~specutils.io.parsing_utils.read_fileobj_or_hdulist`.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
          FITS file name, object (provided from name by Astropy I/O Registry),
          or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    data: SpectrumList
        The spectrum that is represented by the 'loglam' (wavelength) and 'flux'
        data columns in the BINTABLE extension of the FITS `file_obj`.
    """
    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        header = hdulist[0].header
        meta = {'header': header}

        bunit = header.get('BUNIT', '1e-17 erg / (Angstrom cm2 s)')
        if 'Ang' in bunit and 'strom' not in bunit:
            bunit = bunit.replace('Ang', 'Angstrom')
        flux_unit = Unit(bunit)

        # spectrum is in HDU 1
        flux = hdulist[1].data['flux'] * flux_unit

        uncertainty = InverseVariance(hdulist[1].data['ivar'] / flux_unit**2)

        dispersion = 10**hdulist[1].data['loglam']
        dispersion_unit = Unit('Angstrom')

        mask = hdulist[1].data['and_mask'] != 0

    return Spectrum1D(flux=flux, spectral_axis=dispersion * dispersion_unit,
                      uncertainty=uncertainty, meta=meta, mask=mask)


@data_loader(
    label="DESI coadd", identifier=coadd_identify, extensions=['fits'],
    priority=10,
)
def coadd_loader(file_obj, **kwargs):
    """
    Loader for DESI coadd_ files.

    .. _coadd: https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/coadd-SURVEY-PROGRAM-PIXNUM.html

    Any keyword arguments are passed to
    :func:`~specutils.io.parsing_utils.read_fileobj_or_hdulist`.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
          FITS file name, object (provided from name by Astropy I/O Registry),
          or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    data: SpectrumList
        The spectrum that is represented by the 'loglam' (wavelength) and 'flux'
        data columns in the BINTABLE extension of the FITS `file_obj`.
    """
    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        header = hdulist[0].header
        meta = {'header': header}

        bunit = header.get('BUNIT', '1e-17 erg / (Angstrom cm2 s)')
        if 'Ang' in bunit and 'strom' not in bunit:
            bunit = bunit.replace('Ang', 'Angstrom')
        flux_unit = Unit(bunit)

        # spectrum is in HDU 1
        flux = hdulist[1].data['flux'] * flux_unit

        uncertainty = InverseVariance(hdulist[1].data['ivar'] / flux_unit**2)

        dispersion = 10**hdulist[1].data['loglam']
        dispersion_unit = Unit('Angstrom')

        mask = hdulist[1].data['and_mask'] != 0

    return Spectrum1D(flux=flux, spectral_axis=dispersion * dispersion_unit,
                      uncertainty=uncertainty, meta=meta, mask=mask)


def _read_desi(file_obj, **kwargs):
    """
    Read DESI data from a FITS file.

    This contains the common, low-level code for reading spectra and coadd files.

    Any keyword arguments are passed to
    :func:`~specutils.io.parsing_utils.read_fileobj_or_hdulist`.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
          FITS file name, object (provided from name by Astropy I/O Registry),
          or HDUList (as resulting from astropy.io.fits.open()).
    single : bool
    skip_hdus?
    select_columns?

    Returns
    -------
    SpectrumList
        The data.
    """
    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:

    ftype = np.float64
    # if single:
    #     ftype = np.float32

    infile = os.path.abspath(infile)
    if not os.path.isfile(infile):
        raise IOError("{} is not a file".format(infile))

    t0 = time.time()
    hdus = fitsio.FITS(infile, mode="r")
    nhdu = len(hdus)

    if targetids is not None and rows is not None:
        raise ValueError('Set rows or targetids but not both')

    #- default skip_hdus empty set -> include everything, without
    #- having to check for None before checking if X is in skip_hdus
    if skip_hdus is None:
        skip_hdus = set()

    #- Map targets -> rows and exp_rows.
    #- Note: coadds can have uncoadded EXP_FIBERMAP HDU with more rows than
    #- the coadded FIBERMAP HDU, so track rows vs. exp_rows separately
    exp_rows = None
    if targetids is not None:
        targetids = np.atleast_1d(targetids)
        file_targetids = hdus["FIBERMAP"].read(columns="TARGETID")
        rows = np.where(np.isin(file_targetids, targetids))[0]
        if 'EXP_FIBERMAP' in hdus and 'EXP_FIBERMAP' not in skip_hdus:
            exp_targetids = hdus["EXP_FIBERMAP"].read(columns="TARGETID")
            exp_rows = np.where(np.isin(exp_targetids, targetids))[0]
        if len(rows) == 0:
            return Spectra()
    elif rows is not None:
        rows = np.asarray(rows)
        # figure out exp_rows
        file_targetids = hdus["FIBERMAP"].read(rows=rows, columns="TARGETID")
        if 'EXP_FIBERMAP' in hdus and 'EXP_FIBERMAP' not in skip_hdus:
            exp_targetids = hdus["EXP_FIBERMAP"].read(columns="TARGETID")
            exp_rows = np.where(np.isin(exp_targetids, file_targetids))[0]

    if select_columns is None:
        select_columns = dict()

    for extname in ("FIBERMAP", "EXP_FIBERMAP", "SCORES", "EXTRA_CATALOG"):
        if extname not in select_columns:
            select_columns[extname] = None

    # load the metadata.
    meta = dict(hdus[0].read_header())

    # initialize data objects

    bands = []
    fmap = None
    expfmap = None
    wave = None
    flux = None
    ivar = None
    mask = None
    res = None
    extra = None
    extra_catalog = None
    scores = None

    # For efficiency, go through the HDUs in disk-order.  Use the
    # extension name to determine where to put the data.  We don't
    # explicitly copy the data, since that will be done when constructing
    # the Spectra object.

    for h in range(1, nhdu):
        name = hdus[h].read_header()["EXTNAME"]
        log.debug('Reading %s', name)
        if name == "FIBERMAP":
            if name not in skip_hdus:
                fmap = encode_table(
                    Table(
                        hdus[h].read(rows=rows, columns=select_columns["FIBERMAP"]),
                        copy=True,
                    ).as_array()
                )
        elif name == "EXP_FIBERMAP":
            if name not in skip_hdus:
                expfmap = encode_table(
                    Table(
                        hdus[h].read(rows=exp_rows, columns=select_columns["EXP_FIBERMAP"]),
                        copy=True,
                    ).as_array()
                )
        elif name == "SCORES":
            if name not in skip_hdus:
                scores = encode_table(
                    Table(
                        hdus[h].read(rows=rows, columns=select_columns["SCORES"]),
                        copy=True,
                    ).as_array()
                )
        elif name == "EXTRA_CATALOG":
            if name not in skip_hdus:
                extra_catalog = encode_table(
                    Table(
                        hdus[h].read(
                            rows=rows, columns=select_columns["EXTRA_CATALOG"]
                        ),
                        copy=True,
                    ).as_array()
                )
        else:
            # Find the band based on the name
            mat = re.match(r"(.*)_(.*)", name)
            if mat is None:
                raise RuntimeError(
                    "FITS extension name {} does not contain the band".format(name)
                )
            band = mat.group(1).lower()
            type = mat.group(2)
            if band not in bands:
                bands.append(band)
            if type == "WAVELENGTH":
                if wave is None:
                    wave = {}
                # - Note: keep original float64 resolution for wavelength
                wave[band] = native_endian(hdus[h].read())
            elif type == "FLUX":
                if flux is None:
                    flux = {}
                flux[band] = _read_image(hdus, h, ftype, rows=rows)
            elif type == "IVAR":
                if ivar is None:
                    ivar = {}
                ivar[band] = _read_image(hdus, h, ftype, rows=rows)
            elif type == "MASK" and type not in skip_hdus:
                if mask is None:
                    mask = {}
                mask[band] = _read_image(hdus, h, np.uint32, rows=rows)
            elif type == "RESOLUTION" and type not in skip_hdus:
                if res is None:
                    res = {}
                res[band] = _read_image(hdus, h, ftype, rows=rows)
            elif type != "MASK" and type != "RESOLUTION" and type not in skip_hdus:
                # this must be an "extra" HDU
                log.debug('Reading extra HDU %s', name)
                if extra is None:
                    extra = {}
                if band not in extra:
                    extra[band] = {}

                extra[band][type] = _read_image(hdus, h, ftype, rows=rows)

    hdus.close()
    duration = time.time() - t0
    log.info(iotime.format("read", infile, duration))

    # Construct the Spectra object from the data.  If there are any
    # inconsistencies in the sizes of the arrays read from the file,
    # they will be caught by the constructor.

    spec = Spectra(
        bands,
        wave,
        flux,
        ivar,
        mask=mask,
        resolution_data=res,
        fibermap=fmap,
        exp_fibermap=expfmap,
        meta=meta,
        extra=extra,
        extra_catalog=extra_catalog,
        single=single,
        scores=scores,
    )

    return spec
