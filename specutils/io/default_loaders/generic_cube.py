# ~/.specutils/my_custom_loader.py
#
#  Load a FITS cube , extract the spectrum at the reference pixel, but this
#  can be optionally overriden.
#
#  21-apr-2016  Peter Teuben    hackday at "SPECTROSCOPY TOOLS IN PYTHON WORKSHOP" STSCI

import numpy as np
from astropy.units import Unit
from astropy.wcs import WCS

from ...spectra import Spectrum1D
from ..parsing_utils import read_fileobj_or_hdulist


# Define an optional identifier. If made specific enough, this circumvents the
# need to add `format="my-format"` in the `Spectrum1D.read` call.
def identify_generic_fits(origin, *args, **kwargs):
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        return (hdulist[0].header['NAXIS'] == 3)


# not yet ready because it's not generic enough and does not use column_mapping
# @data_loader("Cube", identifier=identify_generic_fits, extensions=['fits'])
def generic_fits(file_obj, **kwargs):

    with read_fileobj_or_hdulist(file_obj, **kwargs) as hdulist:
        header = hdulist[0].header
        data3  = hdulist[0].data
        wcs    = WCS(header)
        shape  = data3.shape

        # take the reference pixel if the pos= was not supplied by the reader
        if 'pos' in kwargs:
            ix = kwargs['pos'][0]
            iy = kwargs['pos'][1]
        else:
            ix = int(wcs.wcs.crpix[0])
            iy = int(wcs.wcs.crpix[1])

        # grab a spectrum from the cube
        if len(shape) == 3:
            data = data3[:,iy,ix]
        elif len(shape) == 4:
            data = data3[:,:,iy,ix].squeeze()
            # make sure this is a 1D array
            # if len(data.shape) != 1:
            #    raise Exception("not a true cube")
        else:
            raise ValueError(f"Unexpected shape {shape}.")

        # store some meta data
        meta = {'header': header}
        meta['xpos'] = ix
        meta['ypos'] = iy

        # attach units (get it from header['BUNIT'] - what about 'JY/BEAM '
        #  NOTE:  astropy doesn't support beam, but see comments in radio_beam
        data = data * Unit("Jy")

        # now figure out the frequency axis....
        sp_axis = 3
        naxis3 = header['NAXIS%d' % sp_axis]
        cunit3 = wcs.wcs.cunit[sp_axis-1]
        crval3 = wcs.wcs.crval[sp_axis-1]
        cdelt3 = wcs.wcs.cdelt[sp_axis-1]
        crpix3 = wcs.wcs.crpix[sp_axis-1]

        freqs = np.arange(naxis3) + 1
        freqs = (freqs - crpix3) * cdelt3 + crval3

        freqs = freqs * cunit3

        # should wcs be transformed to a 1D case ?

    return Spectrum1D(flux=data, wcs=wcs, meta=meta, spectral_axis=freqs)
    # return Spectrum1D(flux=data, wcs=wcs, meta=meta)     # this does not work yet
