import os

from specutils import Spectrum1D



def ascii_reader(filename, filter, **kwargs):
    """Like :func:`fits_reader` but for ASCII file."""
    name = os.path.basename(filename.name.rstrip(os.sep)).rsplit('.', 1)[0]
    tab = ascii.read(filename, **kwargs)
    cols = tab.colnames

    return Spectrum1D(name=str(name), data=data, dispersion=dispersion,
                      uncertainty=uncertainty, mask=mask, wcs=wcs,
                      unit=unit, dispersion_unit=disp_unit)


def ascii_identify(origin, *args, **kwargs):
    """Check whether given filename is ASCII.
    This is used for Astropy I/O Registry.

    """
    return (isinstance(args[0], str) and
            args[0].lower().split('.')[-1] in ['txt', 'dat'])