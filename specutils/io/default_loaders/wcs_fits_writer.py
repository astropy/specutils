from __future__ import absolute_import, division

from astropy.table import Table
from ..registers import custom_writer


@custom_writer("wcs-fits")
def tabular_fits(spectrum, file_name, **kwargs):
    flux = spectrum.flux.value
    disp = spectrum.dispersion.value
    meta = spectrum.meta

    tab = Table([disp, flux], names=("dispersion", "flux"), meta=meta)

    tab.write(file_name, format="fits")
