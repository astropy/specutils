import numpy as np

from astropy.io import fits
from astropy.nddata.nduncertainty import StdDevUncertainty

from ..spectrum1d import Spectrum1D

class ApertureException(Exception):

    def __init__(self, apertures):
        self.apertures = apertures
    
    def __str__(self):

        message = "There is more than one spectrum in the requested file. " + \
            "Each spectrum corresponds to one aperture. " + \
            "Please specify the aperture desired with the aperture " + \
            "= argument. The available apertures are:\n\n" + str(list(self.apertures))

        return message

def _find_col_index(ColDefs, colname):
    '''find index for a column with certain name

    This function can be removed or made easier, once astropy.io.fits returns a Table object'''
    ind = (np.array(ColDefs.names) == colname)
    return np.where(ind)[0][0]

def _get_unit(unitstring):
    '''The only unit I have seen is ERG/CM2/S/A and that does not parse
    with astropy. For now, just translate that.'''
    if unitstring == 'ERG/CM2/S/A':
        return 'erg / (cm2 s Angstrom)'
    else:
        return unitstring

def read_IUE_mxlo(filename, aperture = ''):
    '''
    Read low resolution IUE spectra in fits format

    IUE spectra are available from MAST in two different formats:
    IUESIPS and NEWSIPS. While the former is an (now obsolete) intrinsic IUE
    format that requires special software to read, NEWSPIS is a fits file.
    This procedure implements a reader for NEWSIPS files.

    Parameters
    ----------
    filename : string
    aperture : string
        Some IUE files contain spectra taken with different apertures.
        This keyword allows the user to select one of them (an exception
        is raised if multiple apterures are present, but this keyword
        is not set appropriatly).

    Returns
    -------
    spec : Spectrum1D

    http://archive.stsci.edu/iue/mdr_help.html?print=1
    
    Sometimes there are two or more apertures in one file e.g. swp02283.mxlo .
    '''
    hdus = fits.open(filename)
    try:
        tab = hdus[1].data
        meta = hdus[0].header
        wave = tab.WAVELENGTH[0]+np.arange(tab.NPOINTS[0])*tab.DELTAW[0]
        if len(tab) >1:
            if aperture in tab.APERTURE:
                index = np.where(tab.APERTURE == aperture)
                flux = tab.FLUX[index].flatten()
                sigma = tab.SIGMA[index].flatten()
                flags = tab.QUALITY[index].flatten()
            else:
                raise ApertureException(tab.APERTURE) 
        else:
            flux = tab.FLUX.flatten()
            sigma = tab.SIGMA.flatten()
            flags = tab.QUALITY.flatten()
    
    finally:
        hdus.close()

    dispersion_unit =  tab.columns[_find_col_index(tab.columns, 'WAVELENGTH')].unit.lower()
    flux_unit = _get_unit(tab.columns[_find_col_index(tab.columns, 'FLUX')].unit)

    spec = Spectrum1D.from_array(wave, flux,
                                 uncertainty=StdDevUncertainty(sigma),
                                 meta=meta, dispersion_unit=dispersion_unit,
                                 unit=flux_unit, mask=(flags!=0),
                                 flags=flags)
    return spec

