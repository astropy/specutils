import numpy as np

from astropy.io import fits
from astropy.nddata.nduncertainty import StdDevUncertainty
from specutils.spectrum1d import Spectrum1D

def rdech(filename):
    '''This code reads some form of the ech files, that REDUCE writes.
     
    ech files are in fits standard, but instead of a wavelength column in the
    table the fit coefficient form the wavelength solution are saved. 
    
    Translated from IDL REDUCE - simplified by a lot
    see: http://www.astro.uu.se/~piskunov/RESEARCH/REDUCE/index.html

    The IDL version offers a lot of specific options, 
    none of them is implemented here.
    Some of them are not really necessary: If you have only a few orders
    in your echelle file (the most common case), then it does not hurt to read
    the entire file every time. If you only want one order, you can crop it
    later in python.

    Make more consistent for astropy later. Now it is only important that it works.
    
    Parameters
    ----------
    filename: string or anything pyfits.open takes
    
    Returns
    -------
    spec: Spectrum1D
 
        
    To-Do
    -----
    There are some things to add to make this compatible with all ech
    files and all version.
    - read versions before 2.2
    - use GAIN, ZAPIND, ZAPVAL header keywords if given
    - use RADVEL keyword if given or defined as input
    - output units
    - output header as meta
    '''
    hdus = fits.open(filename)
    data = hdus[1].data
    
    wframe = 'none'
    units = 'adu'
    
    if 'WAVE' in data.names:
        wave = mkwave(data['WAVE'][0,:])
        wframe = 'obs'
        if 'BARYCORR' in hdus[0].header:
            corr = hdus[0].header['BARYCORR']
            wave = wave * (1. + corr/ 2.9979246e5)
            print 'performing barycentric correction'
            wframe = 'bary'
            if 'RADVEL' in hdus[0].header:
                radvel = hdus[0].header['RADVEL']
                wave = wave * (1. - radvel/ 2.9979246e5)
                print 'Correcting for radial velocity'
                wframe = 'source'
    if 'CONT' in data.names:
        spec = data['SPEC'] / data['CONT'] #normalize continuum
        units = 'norm'
    else:
        spec = data['SPEC']
    spec = np.reshape(spec[0,:,:].T,(-1,spec.shape[1]))
    
    if 'SIG' in data.names:                #true: have uncertainties
        if 'CONT' in data.names:  
            sig = data['SIG'] / data['CONT']   #renormalize uncertainties
        else:
            sig = data['SIG']
        sig = np.reshape(sig[0,:,:].T,(-1,sig.shape[1]))
    else: 
        sig = None
    
    return Spectrum1D(spec, dispersion = wave, uncertainty = StdDevUncertainty(sig))

        
def mkwave(wvc):
    '''Make wavelength array from fit values in ech file header
    
    Parameters
    ----------
    wvc: np.array
        `WAVE` section of `ech` file data
        
    Returns
    -------
    wave: np.array
        array of wavelength in Ang
    '''
    version = wvc[0]
    ncol = wvc[1]
    nord = wvc[2]
    obase = wvc[3]
    ncross = int(wvc[7])
    coldeg = int(wvc[8])
    orddeg = int(wvc[9])
    coeff = wvc[10:]

    if version < 2.:
        raise NotImplementedError
    elif version < 2.3:
        ic = np.arange(0, ncol*0.01, 0.01)
    else:
        ic = np.arange(0, ncol*0.001, 0.001)

    #Loop over orders, building wavelength array.
    w = np.zeros((ncol,nord), dtype = np.float)   #init wavelength array
    for i in np.arange(nord):	#loop over build orders
        order = (obase + i) / 100. #current order
        nlc = np.polyval(np.hstack((coeff[coldeg:0:-1],[0])), ic)
        nlo = np.polyval(np.hstack((coeff[coldeg+orddeg:coldeg:-1],[0])), order)

        if version < 2.2:  #no cross-terms
            nlx = 0.0
        else:  #cross-terms added
            inx = coldeg + orddeg +1 #start of cross-terms
            if ncross==4:
                nlx = order*ic *( coeff[inx] + ic*coeff[inx+1] + order*coeff[inx+2] + ic*order*coeff[inx+3] )
            elif ncross==6: 
                nlx = order*ic *( coeff[inx] + ic*coeff[inx+1] + order*coeff[inx+2] + ic*order*coeff[inx+3] + ic**2*coeff[inx+4] + order**2*coeff[inx+5])
            else: 
                raise ValueError('mkwave: unexpected number of cross-terms')

        nl = coeff[0] + nlc + nlo + nlx		#build n*lambda
        w[:,i] = nl / np.float(obase + i)	#recover wavelengths
    return w
