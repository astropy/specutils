from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .registries import loader_registry
from .registries import io_registry

import logging
from ..core.data import Data

import os
import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table
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

    # Let Astropy parse the table for us
    tab = Table.read(hdulist[ref.data['hdu']], format='fits')
    cols = tab.colnames

    # Read in wavelength column if it exists
    dispersion = None
    disp_unit = None
    if hasattr(ref, 'dispersion') and ref.dispersion.get('hdu') is not None:
        if ref.dispersion['hdu'] == ref.data['hdu']:
            wavtab = tab
            wavcols = cols
        else:
            wavtab = Table.read(hdulist[ref.dispersion['hdu']], format='fits')
            wavcols = wavtab.colnames

        dispersion = wavtab[wavcols[ref.dispersion['col']]]
        disp_unit = dispersion.unit
        dispersion = dispersion.data
        if hasattr(dispersion, 'mask'):
            dispersion = dispersion.data

    # Read flux column
    data = tab[cols[ref.data['col']]].data
    unit = tab[cols[ref.data['col']]].unit
    if unit is None:
        try:
            unit = Unit(meta['header'].get('BUNIT', ""))
        except ValueError as e:
            logging.warning(e)
            unit = Unit("")

    # Read data mask
    if hasattr(data, 'mask'):
        mask = data.mask
        data = data.data
    else:
        mask = np.zeros(data.shape)
        if hasattr(ref, 'mask') and ref.mask.get('hdu') is not None:
            try:
                mask = hdulist[ref.mask['hdu']].data
            except IndexError:
                logging.warning("Mask extension not valid.")

    # Read flux uncertainty
    uncertainty = np.zeros(data.shape)
    uncertainty_type = 'std'
    if hasattr(ref, 'uncertainty') and ref.uncertainty.get('hdu') is not None:
        if ref.uncertainty['hdu'] == ref.data['hdu']:
            errtab = tab
            errcols = cols
        else:
            errtab = Table.read(hdulist[ref.uncertainty['hdu']], format='fits')
            errcols = errtab.colnames

        try:
            uncertainty = errtab[errcols[ref.uncertainty['col']]]
            if (uncertainty.unit is not None and unit != Unit("") and
                    uncertainty.unit != unit):
                uncertainty = uncertainty.to(unit).value # TODO: Make this robust
            else:
                uncertainty = uncertainty.data
                if hasattr(uncertainty, 'mask'):
                    uncertainty = uncertainty.data
        except AttributeError:
            pass

        uncertainty_type = ref.uncertainty.get('type', 'std')

    # This will be dictated by the type of the uncertainty
    uncertainty = StdDevUncertainty(uncertainty)

    hdulist.close()

    return Data(name=name, data=data, unit=unit, uncertainty=uncertainty,
                mask=mask, wcs=wcs, dispersion=dispersion,
                dispersion_unit=disp_unit)


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
