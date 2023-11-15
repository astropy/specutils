# a hack around tabular-fits that supports putting header info in the primary hdu

from astropy.wcs import WCS
from specutils import Spectrum1D
import astropy.units as u
from astropy.io import fits
import numpy as np
from astropy.table import Table

def mef_tabular_fits_writer(spectrum, file_name='my_func_output.fits', hdu=0, update_header=False, **kwargs):

    # no matter what the dimensionality of 'flux', there will always be a primary header and
    # a single BinTableHDU extension. the `hdu` arg controls where `meta.header` is dumped.
    hdulist = [fits.PrimaryHDU(), None]
    pri_header_initial = hdulist[0].header

    if hdu > 1:
        raise ValueError('`hdu`, which controls which extension `meta.header` is written to, must either be 0 or 1')
        
    # we expect anything that should be written out to the file header to be in `meta.header`
    # This should be a `fits.header.Header` object, so convert to if meta.header is a dictionary
    header = spectrum.meta.get('header', fits.header.Header())  # if no meta.header, create empty `fits.header.Header`

    if hdu == 0:
        header.update(pri_header_initial)  # update with initial values created from primary header

    if not isinstance(header, fits.header.Header):
        if isinstance(header, dict):
            header = fits.header.Header(header)
        else:
            raise ValueError('`Spectrum1d.meta.header must be `fits.header.Header` or dictionary.')

    if update_header:
        hdr_types = (str, int, float, complex, bool,
                     np.floating, np.integer, np.complexfloating, np.bool_)
        header.update([keyword for keyword in spectrum.meta.items() if
                       isinstance(keyword[1], hdr_types)])

    # Strip header of FITS reserved keywords
    for keyword in ['NAXIS', 'NAXIS1', 'NAXIS2']:
        header.remove(keyword, ignore_missing=True)

    # Add dispersion array and unit
    wtype = kwargs.pop('wtype', spectrum.spectral_axis.dtype)
    wunit = u.Unit(kwargs.pop('wunit', spectrum.spectral_axis.unit))

    disp = spectrum.spectral_axis.to(wunit, equivalencies=u.spectral())

    # Mapping of spectral_axis types to header TTYPE1 (no "torque/work" types!)
    dispname = str(wunit.physical_type)
    if dispname == "length":
        dispname = "wavelength"
    elif "energy" in dispname:
        dispname = "energy"

    # Add flux array and unit
    ftype = kwargs.pop('ftype', spectrum.flux.dtype)
    funit = u.Unit(kwargs.pop('funit', spectrum.flux.unit))
    flux = spectrum.flux.to(funit, equivalencies=u.spectral_density(disp))

    columns = [disp.astype(wtype), flux.astype(ftype)]
    colnames = [dispname, "flux"]

    # Include uncertainty - units to be inferred from spectrum.flux
    if spectrum.uncertainty is not None:
        try:
            unc = (
                spectrum
                .uncertainty
                .represent_as(StdDevUncertainty)
                .quantity
                .to(funit, equivalencies=u.spectral_density(disp))
            )
            columns.append(unc.astype(ftype))
            colnames.append("uncertainty")
        except RuntimeWarning:
            raise ValueError("Could not convert uncertainty to StdDevUncertainty due"
                             " to divide-by-zero error.")

    # For > 1D data transpose from row-major format
    for c in range(1, len(columns)):
        if columns[c].ndim > 1:
            columns[c] = columns[c].T

    tab = Table(columns, names=colnames)
    hdulist[1] = fits.BinTableHDU(tab)

    # now figure out where meta.header (and if specified, addl meta keys) should go
    if hdu==0:
        print('here')
        hdulist[0].header.update(header)
    else:  # must otherwise be 1, we already checked this
        hdulist[1].header.update(header)

    hdulist = fits.HDUList(hdulist)

    hdulist.writeto(file_name, overwrite=True)