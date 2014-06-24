# http://iraf.net/irafdocs/specwcs.php
import re

import numpy as np

import astropy.io.fits as fits
from astropy.wcs import WCS
from ..spectrum1d import Spectrum1D

class IrafFormatError(Exception):
    '''
    '''
    pass

'''
Missing important things
onedspec format
LTM, CD keywords, CDELT

'''

### implement several possible dispersion functions
def polynomial(func, spec):
    '''Evaluate polynomial dispersion functions
    
    IRAF defines two possible polynomials: Chebyshev and Legendre
    
    Parameters
    ----------
    func : function
        Possible functions are Chebyshev and Legendre
    spec : list of strings
        This is split from the "specN= A B C ..." string in a fits header,
        that defines the WCS to use
    
    Returns
    -------
    pmin, pmax: integer
        min and max coordinate descibed by this dispersion relation
        in physical pixles
    result : np.array
       wavelengths calculated form the dispersion relation
    '''
    order = int(spec.pop(0))
    pmin = int(float(spec.pop(0)))
    pmax = int(float(spec.pop(0)))
    n = np.linspace(-1,1,pmax-pmin+1)
    temp = np.zeros((pmax-pmin+1, order), dtype = np.float)
    temp[:,0] = 1.
    temp[:,1] = n
    for i in np.arange(2,order):
        temp[:,i] = func(n, i, temp)
    for i in np.arange(0,order):
        temp[:,i] = float(spec.pop(0)) * temp[:,i]
    return pmin, pmax, temp.sum(axis=1)

def chebyshev(n, i, temp):
    '''Iterative definition of Chebyshev polynomial
    
    Implented here from the iraf.noao.onddspec documentation
    see http://iraf.net/irafdocs/specwcs.php
    
    Parameters
    ----------
    n : np.array
        input normalized coordinate array from -1..1
    i : float or integer
        degree of current iteration
    temp : np.array of dimension (n,i)
        temporary array holding previous iterations
    '''
    return 2. * n * temp[:,i-1] - temp[:,i-2]

def legendre(n, i, temp):
    '''Iterative definition of Legendre polynomial
    
    Implented here from the iraf.noao.onddspec documentation
    see http://iraf.net/irafdocs/specwcs.php
    
    Parameters
    ----------
    n : np.array
        input normalized coordinate array from -1..1
    i : float or integer
        degree of current iteration
    temp : np.array of dimension (n,i)
        temporary array holding previous iterations
    '''
    return ((2*i-1)*n*temp[:,i-1]-(i-1)*temp[:,i-2])/i

def linear_spline(spec):
    '''Evaluate linear spline dispersion function
    
    Parameters
    ----------
    spec : list of strings
        This is split from the "specN= A B C ..." string in a fits header,
        that defines the WCS to use
    
    Returns
    -------
    pmin, pmax: integer
        min and max coordinate descibed by this dispersion relation
        in physical pixles
    result : np.array
       wavelengths calculated fomr the dispersion relation
    '''
    npieces = int(spec.pop(0))
    pmin = int(float(spec.pop(0)))
    pmax = int(float(spec.pop(0)))
    s = np.linspace(0, npieces, pmax-pmin+1)
    temp = np.zeros((pmax-pmin+1), dtype = np.float)
    for i in np.arange(0,npieces+1):
        weight = np.clip(1.-np.abs(i-s), 0., 1.)
        temp += weight * float(spec.pop(0))
    return pmin, pmax, temp

def cubic_spline(spec):
    '''Evaluate linear spline dispersion function
    
    Parameters
    ----------
    spec : list of strings
        This is split from the "specN= A B C ..." string in a fits header,
        that defines the WCS to use
    
    Returns
    -------
    pmin, pmax: integer
        min and max coordinate descibed by this dispersion relation
        in physical pixles
    result : np.array
       wavelengths calculated fomr the dispersion relation
    '''
    npieces = int(spec.pop(0))
    pmin = int(float(spec.pop(0)))
    pmax = int(float(spec.pop(0)))
    s = np.linspace(0, npieces, pmax-pmin+1)
    temp = np.zeros((pmax-pmin+1), dtype = np.float)
    for i in np.arange(0,npieces+3):
        c_i = float(spec.pop(0))
        a = np.linspace(0, 1, (pmax-pmin+1)/npieces)
        b = 1. - a
        ind = (s>=i-3) & (s<=i-2)
        if np.sum(ind) > 0: temp[ind] += c_i * b**3
        ind = (s>=i-2) & (s<=i-1)
        if np.sum(ind) > 0: temp[ind] += c_i * (1.+3*b*(1+a*b))
        ind = (s>=i-1) & (s<=i)
        if np.sum(ind) > 0: temp[ind] += c_i * (1.+3*a*(1+a*b))
        ind = (s>=i) & (s<=i+1)
        if np.sum(ind) > 0: temp[ind] += c_i * a**3
    return pmin, pmax, temp



def wave_multispec(header):
    '''calculate wavelength array from MULTISPEC header
    
    This function extracts information from the header of a fits file that
    contains spectral data. From this information, the wavelength array is 
    calculated.
    
    ..note: While some checking on the fits header keywords is done, the
        implementation is not complete. Also, many telescopes use non-standard
        keywords or do not write required values. Be wary!
    
    Parameters
    ----------
    header : astropy.io.fits.header.Header
        header of fits file
    
    Returns
    -------
    wave : np.array
        wavelength array in the format specified in the header. Depending on the 
        header keywords, this can be one or two dimensional.
    '''
    wave = np.empty((header['NAXIS2'], header['NAXIS1']), dtype= np.float)
    wave[:] = np.nan 
    head = header['WAT*']
    formstring = []
    for i in range(len(head)):
        formstring.append(head[i])
    formstring = "".join(formstring)
    formstring = formstring.replace('= ', '=')
    formstring = formstring.replace(' =', '=')
    specs = re.findall(r'spec[0-9]+="[0-9eE .+-]*"', formstring)
    for spec in specs:
        spec=spec[4:].replace('=',' ').replace('"','')
        #specN = ap beam dtype w1 dw nw z aplow aphigh [functions_i]
        #    0    1  2     3    4  5  6 7   8     9     10:
        spec = spec.split()
        N = int(spec.pop(0))
        ap = spec.pop(0)
        beam = spec.pop(0)
        dtype = int(spec.pop(0))
        w1 = float(spec.pop(0))
        dw = float(spec.pop(0))
        nw = int(spec.pop(0))
        z = float(spec.pop(0))
        aplow = float(spec.pop(0))
        aphigh = float(spec.pop(0))
        if dtype == -1:   # no dispersion solution in header
            pass
        elif dtype == 0:  # linear dispersion
            stop = w1 + dw * (nw +0.5)
            wave[N-1,:] = np.arange(w1, stop, dw)/(1.+z)
        elif dtype == 1:  # log-linear dispersion
            stop = w1 + dw * (nw +0.5)
            wave[N-1,:] = 10.**np.arange(w1, stop, dw)/(1.+z)
        elif dtype == 2:  # non-linear dispersion
            wave[N-1,:] = 0
            while len(spec) > 0:
                weight = float(spec.pop(0))
                offset = float(spec.pop(0))
                func = int(spec.pop(0))
                if func == 1:
                    pmin, pmax, thiswave = polynomial(chebyshev, spec)
                elif func == 2:
                    pmin, pmax, thiswave = polynomial(legendre, spec)
                elif func == 3:
                    pmin, pmax, thiswave = cubic_spline(spec)
                elif func == 4:
                    pmin, pmax, thiswave = linear_spline(spec)
                elif func == 5:
                    # TBD: pixel coordinate array
                    raise NotImplementedError
                elif func == 6:
                    #TBD: sampled coordinate array
                    raise NotImplementedError
                else:
                    raise IrafFormatError('File does not conform to IRAF format - non-linear of dispersion function must be [1,2,3,4,5,6]')
                wave[N-1,pmin-1:pmax] += weight * (offset + thiswave)
            wave[N-1,:] = wave[N-1,:]/(1.+z)

        else:
            raise IrafFormatError('File does not conform to IRAF format - type of dispersion function must be [-1,0,1,2]')
    return wave

def wave_1dspec(header):
    '''calculate wavelength array from header for 1 dim WCS
    
    This function extracts information from the header of a fits file that
    contains spectral data. From this information, the wavelength array is 
    calculated.
    
    Parameters
    ----------
    header : astropy.io.fits.header.Header
        header of fits file
    
    Returns
    -------
    wave : np.array
        wavelength array in the format specified in the header.
    '''
    w = WCS(header)
    return w.wcs_pix2world(np.arange(header['NAXIS1'])[:,np.newaxis], 0).squeeze()


def make_IRAF_wave(header):
    old_STRIP_HEADER_WHITESPACE = fits.STRIP_HEADER_WHITESPACE()
    try:
        fits.STRIP_HEADER_WHITESPACE.set(False)
        if header['NAXIS'] ==2:
            if header['WCSDIM'] == 2:
                if header['CTYPE1'] == 'MULTISPE' or header['CTYPE2'] == 'MULTISPE':
                    if header['CTYPE1'] == 'MULTISPE' and header['CTYPE2'] == 'MULTISPE':
                        wave = wave_multispec(header)
                    else:
                        raise IrafFormatError('File does not conform to IRAF format - multispec spectra must have CTYPE1 and CTYPE2 = MULTISPE')
            else:
                raise NotImplementedError
        elif header['NAXIS'] ==1:
            wave = wave_1dspec(header)
        else:
            raise NotImplementedError
        return wave
    finally:
        fits.STRIP_HEADER_WHITESPACE.set(old_STRIP_HEADER_WHITESPACE)

def read_IRAF_spec(filename):
    '''
    Read spectrum in IRAF WCS=MULTISPEC form
  
    Parameters
    ----------
    filename : string
  
    Returns
    -------
    spec : Spectrum1D

    To-Do
    -----
    - take care of units
    - look for uncertainties
    '''

    hdus = fits.open(filename)
    try:
        flux = hdus[0].data
        meta = hdus[0].header
        wave = make_IRAF_wave(meta)
    
    finally:
        hdus.close()

    spec = Spectrum1D.from_array(wave, flux, meta = meta)
    return spec




