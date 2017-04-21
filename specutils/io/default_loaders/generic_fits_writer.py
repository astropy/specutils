from astropy.table import Table
from specutils.io.registers import custom_writer


@custom_writer("generic-fits-writer")
def generic_fits(spectrum, file_name, **kwargs):
    flux = spectrum.flux.value
    disp = spectrum.dispersion.value
    meta = spectrum.meta

    tab = Table([disp, flux], names=("dispersion", "flux"), meta=meta)

    tab.write(file_name, format="fits")