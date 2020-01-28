"""
Loader for PFS spectrum files.

https://github.com/Subaru-PFS/datamodel/blob/master/datamodel.txt
"""
import os
import re

from astropy.io import fits
from astropy.units import Unit
from astropy.nddata import StdDevUncertainty

import numpy as np

from specutils.io.registers import data_loader
from specutils import Spectrum1D

__all__ = ['spec_identify', 'spec_loader']

# This RE matches the file name pattern defined in Subaru-PFS' datamodel.txt :
# "pfsObject-%05d-%s-%3d-%08x-%02d-0x%08x.fits" % (tract, patch, catId, objId,
#                                                  nVisit % 100, pfsVisitHash)
_spec_pattern = re.compile(r'pfsObject-(?P<tract>\d{5})-(?P<patch>.{3})-'
                           r'(?P<catId>\d{3})-(?P<objId>\d{8})-'
                           r'(?P<nVisit>\d{2})-(?P<pfsVisitHash>0x\w{8})'
                           r'\.fits')


def spec_identify(origin, *args, **kwargs):
    """
    Check whether given file is FITS. This is used for Astropy I/O Registry.
    """
    filepath = args[0]
    fileobj = None
    if isinstance(args[0], str):
        try:
            fileobj = open(filepath, mode='rb')
        except FileNotFoundError:
            fileobj = None
    elif fits.util.isfile(args[0]):
        fileobj = args[0]
        filepath = fileobj.name
    # Check for `urlopen` object - can only probe content if seekable
    elif hasattr(args[0], 'url') and hasattr(args[0], 'seekable'):
        filepath = args[0].url
        if args[0].seekable():
            fileobj = args[0]

    return (_spec_pattern.match(os.path.basename(filepath)) is not None and
            fits.connect.is_fits(origin, filepath, fileobj, *args))


@data_loader(label="Subaru-pfsObject", identifier=spec_identify,
             extensions=['fits'])
def spec_loader(file_obj, **kwargs):
    """
    Loader for PFS combined spectrum files.

    Parameters
    ----------
    file_obj: str or file-like
        FITS file name or object (provided from name by Astropy I/O Registry).

    Returns
    -------
    data: Spectrum1D
        The spectrum that is represented by the data in this table.
    """
    if isinstance(file_obj, str):
        file_name = file_obj
    else:
        file_name = file_obj.name

    m = _spec_pattern.match(os.path.basename(file_name))

    with fits.open(file_obj, **kwargs) as hdulist:
        header = hdulist[0].header
        meta = {'header': header,
                'tract': m['tract'],
                'patch': m['patch'],
                'catId': m['catId'],
                'objId': m['objId'],
                'nVisit': m['nVisit'],
                'pfsVisitHash': m['pfsVisitHash']}

        # spectrum is in HDU 2
        data = hdulist[2].data['flux']
        unit = Unit('nJy')

        error = hdulist[2].data['fluxVariance']
        uncertainty = StdDevUncertainty(np.sqrt(error))

        wave = hdulist[2].data['lambda']
        wave_unit = Unit('nm')

        mask = hdulist[2].data['mask'] != 0

    return Spectrum1D(flux=data * unit,
                      spectral_axis=wave * wave_unit,
                      uncertainty=uncertainty,
                      meta=meta,
                      mask=mask)
