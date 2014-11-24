# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module contains various functions for spectra
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np


# From Neil Crighton's Barak
def convolve_constant_dv(wa, fl, wa_dv=None, npix=4., vfwhm=None, dv_const=False,
                         bound='extend'):
    """ Convolve a wavelength array with a gaussian of constant
    velocity width.

    If `vfwhm` is specified an intermediate wavelength array with
    constant velocity pixel width is calculated. Otherwise, both
    `wa_dv` and `npix` must be given -- this is faster because no
    intermediate array needs to be calculated.

    Parameters
    ----------
    fl, wa : arrays of floats, length N
      The array to be convolved and its wavelengths.
    vfwhm : float, optional
      Full width at half maximum in velocity space (km/s) of the
      gaussian kernel with which to convolve `fl`.
    npix : float, default 4  [FWHM]
      Number of pixels corresponding to `vfwhm` in `wa_dv` if given,
      otherwise `wa` is interpolated to an array with velocity pixel
      width = vfwhm / npix.
    wa_dv : array of floats, default `None`
      Wavelength array with a constant velocity width (this can be
      generated with make_constant_dv_wa_scale()).
    dv_const : Bool, default 'False'
      Specifies whether the input array has constant velocity
    bound : str ('extend')
      Boundary condition to convolve()

    Returns
    -------
    fl_out : array of length N
      fl convolved with the gaussian kernel with the specified FWHM.     

    Examples
    --------
    >>> from barak import sed
    >>> wa = np.arange(5000, 7000, 0.03)
    >>> fl = np.random.randn(len(wa)) + 1.
    >>> instrument_profile_fwhm = 5.5 # resolution fwhm km/s
    >>> flsmooth = convolve_constant_dv(wa, fl, vfwhm=instrument_profile)

    To avoid generating a new constant velocity scale for multiple
    convolutions of the same spectrum, you can create your own and
    pass it as the wa_dv argument:

    >>> from barak import sed
    >>> subpixkms = 1.3  # km/s    
    >>> wa_dv = sed.make_constant_dv_wa_scale(wa[0], wa[-1], subpixkms)
    >>> npix = instrument_profile / subpixkms
    >>> flsmooth = convolve_constant_dv(wa, fl, wa_dv=wa_dv, npix=npix)
    """
    from astropy.convolution import convolve, Gaussian1DKernel
    
    # interpolate to the log-linear scale, convolve, then
    # interpolate back again.
    if vfwhm is not None:
        wa_dv = make_constant_dv_wa_scale(wa[0], wa[-1], float(vfwhm)/npix)
    if dv_const is False:
        fl_dv = np.interp(wa_dv, wa, fl)
    else: fl_dv = fl
    # Kernel
    stddev = npix / (2*np.sqrt(2*np.log(2)))
    kernel = Gaussian1DKernel(stddev=stddev)
    # Convolve
    fl_dv_smoothed = convolve(fl_dv, kernel, boundary=bound)
    # Interpolate as needed
    if dv_const is False:
        fl_out = np.interp(wa, wa_dv, fl_dv_smoothed)
    else: fl_out = fl_dv_smoothed
    # Return
    return fl_out

# From Neil Crighton's Barak
def make_constant_dv_wa_scale(wmin, wmax, dv):
    """ Make a wavelength scale with bin widths corresponding to a
    constant velocity width scale.

    Parameters
    ----------
    wmin, wmax : floats
      Start and end wavelength.
    dv : float
      Velocity pixel width.

    Returns
    -------
    wa : ndarray
      Wavelength scale. 

    See Also
    --------
    barak.spec.make_wa_scale
    """
    dlogw = np.log10(1 + dv/c_kms)
    # find the number of points needed.
    npts = int(np.log10(wmax / wmin) / dlogw)
    wa = wmin * 10**(np.arange(npts)*dlogw)
    return wa
