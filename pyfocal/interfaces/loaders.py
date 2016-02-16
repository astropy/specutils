from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .registries import loader_registry
from .registries import io_registry

import logging
from ..core.data import Data

import os
import numpy as np
from astropy.io import ascii, fits
from astropy.wcs import WCS
from astropy.nddata import StdDevUncertainty
from astropy.units import Unit


def fits_reader(filename, filter, **kwargs):
    """
    This generic function will query the loader factory, which has already
    loaded the yaml configuration files, in an attempt to parse the
    associated fits file.
    """
    logging.info("Attempting to open '{}' using filter '{}'.".format(
            filename, filter))
    #print("Attempting to open '{}' using filter '{}'.".format(
    #        filename, filter))

    name = os.path.basename(filename.name.rstrip(os.sep)).rsplit('.', 1)[0]
    hdulist = fits.open(filename, **kwargs)
    ref = loader_registry.get(filter)

    meta = ref.meta
    header = dict(hdulist[ref.wcs['hdu']].header)
    meta['header'] = header
    wcs = WCS(hdulist[ref.wcs['hdu']].header)
    data = hdulist[ref.data['hdu']].data
    uncertainty = np.zeros(data.shape)
    uncertainty_type = 'std'

    try:
        unit = Unit(meta['header'].get('BUNIT', ""))
    except ValueError as e:
        logging.warning(e)
        unit = Unit("")

    mask = np.zeros(data.shape)

    # Grab the data array, use columns if necessary
    if ref.data.get('col') is not None:
        try:
            data = data[data.columns[ref.data['col']].name]
        except AttributeError:
            logging.warning("No such columns available.")

    if hasattr(ref, 'uncertainty') and ref.uncertainty.get('hdu') is not None:
        try:
            uncertainty = hdulist[ref.uncertainty['hdu']].data

            if ref.uncertainty.get('col') is not None:
                uncertainty = uncertainty[
                    uncertainty.columns[ref.uncertainty['col']].name]
        except AttributeError:
            pass

        uncertainty_type = ref.uncertainty.get('type', 'std')

    # This will be dictated by the type of the uncertainty
    uncertainty = StdDevUncertainty(uncertainty)

    if hasattr(ref, 'mask') and ref.mask.get('hdu') is not None:
        try:
            mask = hdulist[ref.mask['hdu']].data
        except IndexError:
            logging.warning("Mask extension not valid.")

    hdulist.close()

    return Data(name=name, data=data, uncertainty=uncertainty, mask=mask,
                wcs=wcs, unit=unit)


def fits_identify(origin, *args, **kwargs):
    """Check whether given filename is FITS."""
    return (isinstance(args[0], str) and
            args[0].lower().split('.')[-1] in ['fits', 'fit'])


def ascii_reader(filename, filter, **kwargs):
    """Parse given ASCII file to spectrum data."""
    name = os.path.basename(filename.name.rstrip(os.sep)).rsplit('.', 1)[0]
    tab = ascii.read(filename, **kwargs)
    cols = tab.colnames
    ref = loader_registry.get(filter)

    meta = ref.meta
    meta['header'] = {}

    # Only loads KEY=VAL comment entries into header
    if 'comments' in tab.meta:
        for s in tab.meta['comments']:
            if '=' not in s:
                continue
            s2 = s.split('=')
            meta['header'][s2[0]] = s2[1]

    wcs = None
    wave = tab[cols[ref.dispersion['col']]]
    dispersion = wave.data
    flux = tab[cols[ref.data['col']]]
    data = flux.data
    uncertainty = np.zeros(data.shape)
    uncertainty_type = 'std'

    if flux.unit is None:
        unit = Unit(ref.data.get('unit', 'erg / (Angstrom cm2 s)'))
    else:
        unit = flux.unit

    if wave.unit is None:
        disp_unit = Unit(ref.dispersion.get('unit', 'Angstrom'))
    else:
        disp_unit = wave.unit

    mask = np.zeros(data.shape)

    if hasattr(ref, 'uncertainty') and ref.uncertainty.get('col') is not None:
        try:
            uncertainty = tab[cols[ref.uncertainty['col']]].data
        except IndexError:
            pass  # Input has no uncertainty column
        else:
            uncertainty_type = ref.uncertainty.get('type', 'std')

    # This will be dictated by the type of the uncertainty
    uncertainty = StdDevUncertainty(uncertainty)

    if hasattr(ref, 'mask') and ref.mask.get('col') is not None:
        try:
            mask = tab[cols[ref.mask['col']]].data
        except IndexError:
            pass  # Input has no mask column

    return Data(name=name, data=data, dispersion=dispersion,
                uncertainty=uncertainty, mask=mask, wcs=wcs, unit=unit,
                dispersion_unit=disp_unit)


def ascii_identify(origin, *args, **kwargs):
    """Check whether given filename is ASCII."""
    return (isinstance(args[0], str) and
            args[0].lower().split('.')[-1] in ['txt', 'dat'])


# Add IO reader/identifier to io registry
io_registry.register_reader('fits', Data, fits_reader)
io_registry.register_identifier('fits', Data, fits_identify)

io_registry.register_reader('ascii', Data, ascii_reader)
io_registry.register_identifier('ascii', Data, ascii_identify)
