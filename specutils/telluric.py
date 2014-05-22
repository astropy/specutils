from astropy import log
from astropy.table import Table, Column
from astropy.nddata import NDData, FlagCollection
from astropy.io import fits
from astropy.modeling import models
from astropy.utils import misc

from .wcs import Spectrum1DLookupWCS, Spectrum1DLinearWCS
from specutils import spectrum1d, rvmeasure
import numpy as np
import os
from scipy import signal,interpolate

from pylab import *

def fixgstar(wave,flux,template=None,debug=False,order=3,smooth=1.0,fitRange=None,nPixFit=2,saveFile=None):
    '''
    This function takes in a spectrum of a G2V star and divides by the
    solar spectrum to remove the stellar lines to make a telluric
    spectrum to remove atmospheric lines.

    template -- the template to use (by default is data/solarspectrum.hdf5
                (can be .txt instead if passed in)
    order -- the order of the polynomial fit to normalize the input spectrum
    smooth -- additional factor to smooth the template spectrum beyond just
              the convolution with the resolution (default = 1.0). 
    '''

    if template is None:
        template = os.path.join(os.path.dirname(__file__), 'data/solarspectrum.hdf5')

    print 'Using template: '+template
    
    if os.path.isfile(template):
        if os.path.splitext(template)[1] == '.hdf5':
            tab = Table.read(template,path='data')
            solWave = tab['wave']
            solFlux = tab['flux']
            solWave = solWave*10.0 # convert from nm to Angstroms
        else:
            solWave, solFlux = np.loadtxt(template,unpack=True)

        if debug:
            clf()
            ## psf = models.Gaussian1DModel(mean=0.0,stddev=5.0,amplitude=1.0)
            ## psf = psf(np.arange(-20,20))
            ## psf = psf/np.sum(psf)
            ## solFlux = signal.fftconvolve(solFlux,psf,mode='same')
            ## semilogy(solWave,solFlux)

        # cut out the regions of the template that are not within the input range
        goodWaves = np.where((solWave >= wave[0]) & (solWave <= wave[-1]))[0]

        # take out a few points more so that we can interpolate later
        if goodWaves[0] > 5:
            goodWaves = np.append(goodWaves[0]-np.arange(1,5),goodWaves)
        if goodWaves[-1] < len(solWave)+5:
            goodWaves = np.append(goodWaves,goodWaves[-1]+np.arange(1,5))

        solWave = solWave[goodWaves]
        solFlux = solFlux[goodWaves]

        # lets convolve to the right resolution
        delta = np.mean(np.diff(wave))/np.mean(np.diff(solWave))*smooth
        if delta > 1.0:
            if debug:
                print 'delta: '+str(delta)
            psf = models.Gaussian1D(mean=0.0,stddev=delta,amplitude=1.0)
            psf = psf(np.arange(-np.fix(delta*3.0),np.fix(delta*3)))
            psf = psf/np.sum(psf)
            solFlux = signal.fftconvolve(solFlux,psf,mode='same')
            


        goodIn1 = np.where(flux > 0)[0]
        pFit1 = np.polyfit(wave[goodIn1],flux[goodIn1],order)

        # don't fit the beginning and end points in case there are convolution artifacts
        pFit2 = np.polyfit(solWave[1:-2],solFlux[1:-2],order)
        
        if debug:
            clf()
            subplot(211)
            plot(wave,flux/np.polyval(pFit1,wave),label='Orig.')
            plot(solWave,solFlux/np.polyval(pFit2,solWave),label='Ref.')
            legend()

        fluxTemp = flux/np.polyval(pFit1,wave)
        solFluxTemp = solFlux/np.polyval(pFit2,solWave)
        
        # measure the relative shift of the spectra
        shiftVel, shiftPeak, logInt = rvmeasure.rvshift(wave,fluxTemp,solWave,solFluxTemp,debug=False,fitRange=fitRange,nPixFit=nPixFit)

        if debug:
            print 'velocity shift: ',shiftVel

        # shift the solar spectrum. Should use the original spectrum
        # because we need to remove the BB shape from the standard
        # star too.
        
        solFluxShift = rvmeasure.shiftSpec(solWave,solFlux,-shiftVel)/np.median(solFlux)

        # interpolate the reference spectrum
        solInterp = interpolate.interp1d(solWave,solFluxShift)
        solSpec = solInterp(wave)

        # remove the intrinsic stellar lines and put back the
        # polynomial fit that was used to flatten the spectrum
        correctedFlux = flux/solSpec*np.polyval(pFit1,wave)
        correctedFlux = correctedFlux/np.median(correctedFlux)
        
        if debug:
            subplot(212)
            plot(wave,flux/np.median(flux),label='Orig.')
            plot(wave,solSpec,label='Shift Ref.')
            plot(wave,correctedFlux,label='Final')
            legend()


        specObj = spectrum1d.Spectrum1D.from_array(wave,correctedFlux)
        if saveFile is not None:
            specObj.to_fits(saveFile)
            
        return specObj
    else:
        print 'Template not found: '+template

def fixastar(wave,flux,template=None,**kwargs):
    '''
    Has the same keywords as the fixgstar except uses the vega template instead
    '''
    template = os.path.join(os.path.dirname(__file__), 'data/vega_all.hdf5')
    return fixgstar(wave,flux,template=template,**kwargs)

def test_fixastar():
    outputDir = '/u/tdo/mosdrp/reduced/mercer23/'
    tellFile = outputDir+'astar_hip98640.fits'
    outfile = outputDir+'telluric_standard_130902.fits'

    specObj = spectrum1d.Spectrum1D.from_fits(tellFile)
    flux = specObj.flux
    wave = specObj.dispersion
    good = np.where(flux > 0)[0]

    pFit = np.polyfit(wave[good],flux[good],3)
    print pFit
    clf()
    subplot(2,1,1)
    plot(specObj.dispersion,specObj.flux)
    plot(specObj.dispersion,np.polyval(pFit,wave))
    subplot(2,1,2)
    plot(specObj.dispersion,flux/np.polyval(pFit,wave))

    fixastar(wave,flux,debug=True,smooth=1.0,fitRange=[16200,16800])

def test_fixgstar():
    outputDir = '/u/tdo/mosdrp/reduced/mercer23/'
    tellFile = outputDir+'gstar_hip97003.fits'
    outfile = outputDir+'telluric_standard_130902.fits'
    
    specObj = spectrum1d.Spectrum1D.from_fits(tellFile)
    flux = specObj.flux
    wave = specObj.dispersion
    good = np.where(flux > 0)[0]
    
    pFit = np.polyfit(wave[good],flux[good],3)
    print pFit
    clf()
    subplot(2,1,1)
    plot(specObj.dispersion,specObj.flux)
    plot(specObj.dispersion,np.polyval(pFit,wave))
    subplot(2,1,2)
    plot(specObj.dispersion,flux/np.polyval(pFit,wave))

    fixgstar(wave,flux,debug=True,smooth=1.0,fitRange=[15800,16800],saveFile=outfile)
