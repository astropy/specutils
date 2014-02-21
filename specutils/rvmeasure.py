import numpy as np
from scipy import signal
from scipy import interpolate
import pylab as pl

def gaussian(p, x, normalize = False):
    ''' gaussian model
    Must have a minimum of 3 terms (scale, centroid, sigma)
    p[0] -- scale factor
    p[1] -- centroid
    p[2] -- sigma
    p[3] -- offset
    p[4] -- slope
    '''

    # check to see how many parameters there are in the Gaussian
    if len(p) < 5:
        slope = 0.0
    else:
        slope = p[4]

    if len(p) < 4:
        offset = 0.0
        slope = 0.0
    else:
        offset = p[3]

    if normalize is True:
        scaleFactor = 1.0/(p[1]*np.sqrt(2.0*np.pi))
    else:
        scaleFactor = p[0]
    u = (x - p[1])/p[2]
    
    return scaleFactor*np.exp(-0.5*u*u) + offset + slope*x



def rvshift(wave1, spec1, wave2, spec2, r1 = None, r2 = None,debug=False,
            oversample=1,fitRange=None,lagRange=[-20,20],nPixFit = 3):
    """Measure the relative radial velocity of two spectra using cross
    correlation. Will return the shift of the 2nd spectrum relative to
    the first.


    Parameters:
    -----------
    wave1 : wavelength for the first spectrum
    spec1 : the first spectrum
    wave2 : wavelengths for the second spectrum
    spec2 : the second spectrum


    Keywords:
    -----------
    r1 : spectral resolution of the first spectrum (default: None)
    r2 : spectral resolution of the second spectrum (default: None)
         If r1 and r2 are different, then will convolve the higher
         resolution spectrum to that of the lower resolution one.
    nPixFit: number of points around either side of the correlation
             peak to fit (default: 3).
    
    Returns:
    ---------
    Return an array of [velocity, pixel shift]
    """

    # make a copy of the input spectra
    waveIn1 = np.copy(wave1)
    specIn1 = np.copy(spec1)
    waveIn2 = np.copy(wave2)
    specIn2 = np.copy(spec2)
    if debug:
        pl.clf()
        pl.subplot(121)
        pl.plot(waveIn1,specIn1/np.mean(specIn1))
        pl.plot(waveIn2,specIn2/np.mean(specIn2))
    if (r1 is not None) and (r2 is not None):
        if r1 > r2:
            delt = wave1[1]/(r2*(wave1[1] - wave1[0]))/(2*np.sqrt(2*np.log(2)))
            npsfPix = 4*int(delt)+1
            psf = gaussian([1.0,npsfPix/2,delt,0],np.arange(npsfPix),normalize=True)
            if debug:
                print np.shape(psf)
            #psf = psf/np.sum(psf)
            specIn1 = signal.fftconvolve(specIn1,psf,mode='same')
            if debug:
                pl.plot(waveIn1,specIn1/np.mean(specIn1))
        if r2 > r1:
            delt = wave2[1]/(r1*(wave2[1] - wave2[0]))/(2*np.sqrt(2*np.log(2)))
            npsfPix = 4*int(delt)+1
            psf = gaussian([1.0,npsfPix/2,delt,0],np.arange(npsfPix),normalize=True)
            if debug:
                print np.shape(psf)
            #psf = psf/np.sum(psf)
            specIn2 = signal.fftconvolve(specIn2,psf,mode='same')
            if debug:
                pl.plot(waveIn2,specIn2/np.mean(specIn1))

    if (len(waveIn2) != len(waveIn1)) or (np.sum(waveIn2 - waveIn1) != 0):
        # interpolate the wavelength locations of the second spectrum
        # to the location of the first they are not the same
        
        f = interpolate.interp1d(waveIn2,specIn2)
        specIn2 = f(waveIn1)
        waveIn2 = waveIn1
        if debug:
            pl.plot(waveIn2,specIn2/np.mean(specIn2))
                 
    # convert the wavelengths to be linearly sampled in log
    lWave1 = np.log(waveIn1)
    nW = len(lWave1)
    logWave1 = np.linspace(lWave1[0],lWave1[-1],num=nW*oversample)
    logInt = logWave1[1] - logWave1[0]  # size of each log bin
    linWave1 = np.exp(logWave1) # the normal wavelength corrsponding to
                                      # the new log array
    # the wavelength samples of both spectra should be the same now
    if debug:
        print 'wave1 start, end', wave1[0],wave1[-1]
        print 'interp wave: ',linWave1[0],linWave1[-1]
    logInterp1 = interpolate.interp1d(waveIn1,specIn1,bounds_error=False)
    logInterp2 = interpolate.interp1d(waveIn1,specIn2,bounds_error=False)
    logSpec1 = logInterp1(linWave1)
    logSpec2 = logInterp2(linWave1)

    # choose only a portion of the wavelength region to cross correlate
    if fitRange is not None:
        good = np.where((linWave1 >= fitRange[0]) & (linWave1 <= fitRange[1]))[0]
    else:
        good = np.arange(len(linWave1))
    
    corr = signal.correlate(logSpec1[good]-np.median(logSpec1[good]),logSpec2[good]-np.median(logSpec2[good]),mode='same')
    lags = np.arange(len(good))-len(good)/2

    if debug:
        pl.subplot(122)
        pl.plot(lags,corr)


    # peak velocity corresponding to the pixel peak
    peakInd = np.argmax(corr)
    peakPixVel = -(np.exp(lags[peakInd]*logInt)-1)*3e5

    # fit a polynomial near the peak of the correlation

    pFit = np.polyfit(lags[peakInd-nPixFit:peakInd+nPixFit+1],corr[peakInd-nPixFit:peakInd+nPixFit+1],2)
    shiftPeak = -pFit[1]/(2*pFit[0])
    
    shiftPeakVel = -(np.exp(shiftPeak*logInt)-1)*3e5
    
    if debug:
        print 'logInt: '+str(logInt)
        print 'shift peak: '+str(lags[peakInd])
        print 'shift peak pixel vel: '+str(peakPixVel)
        print 'polyfit: ',pFit
        print 'shift peak fitted: '+str(shiftPeak)
        print 'shift peak fitted vel: '+str(shiftPeakVel)
        pl.plot([lags[peakInd],lags[peakInd]],[0,corr[peakInd]])        
        pl.xlim(lagRange[0],lagRange[1])
        pl.xlabel('Lag (pixels)')
        polyFun = np.poly1d(pFit)
        pl.plot(lags[peakInd-nPixFit:peakInd+nPixFit+1],polyFun(lags[peakInd-nPixFit:peakInd+nPixFit+1]),'r')
    return np.array([shiftPeakVel,shiftPeak,logInt])

def shiftSpec(wave,flux,vel):
    # takes a the wavelengths and spectra and shift by some velocity value. Then resample back into the same grid.
    # wave -- in angstroms
    # velocity -- in km/s

    # shift the velocity
    actualWave = wave*(vel/3e5+1.0)

    # linearly interpolate the spectrum back to the original wavelength locations
    interp = interpolate.interp1d(actualWave,flux,bounds_error=False)
    
    return interp(wave)
    
def test_rvmeasure():
    wave1 = np.linspace(2.0,2.45,num=2000)
    delt = wave1[1]-wave1[0]
    lines1 = np.array([2.1661, 2.20,2.3])
    testRV = 301.0
    lines2 = testRV/3e5*lines1 + lines1
    spec1 = np.zeros(len(wave1))
    spec2 = np.zeros(len(wave1))
    for i in np.arange(len(lines1)):
        spec1 = spec1+gaussian([1.0,lines1[i],4.0*delt,0],wave1)
        spec2 = spec2+gaussian([1.0,lines2[i],2.0*delt,0],wave1)

    rv = rvshift(wave1,spec1+1.0,wave1,spec2+1.0,debug=True,r1=2000,r2=4000,lagRange=[-100,100],fitRange=[2.1,2.4])
    shifted = shiftSpec(wave1*1e4,spec2,-rv[0])
    pl.clf()
    pl.plot(wave1,spec1,label='Ref')
    pl.plot(wave1,spec2,label='Test Spec')
    pl.plot(wave1,shifted,label='Shifted Fit')
    pl.legend()
