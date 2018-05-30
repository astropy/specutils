import astropy.units as u
from astropy.modeling import models, fitting
import numpy as np

from specutils.tests.spectral_examples import simulated_spectra
from specutils.spectra.spectrum1d import Spectrum1D
from specutils.fitting.fitmodels import fit_models, fit_models_simple

def test_fitting_simple_slsqplsq(simulated_spectra):
    """
    This test fits the the first simulated spectrum from the fixture.  The
    initial guesses are manually set here with bounds that essentially make
    sense as the functionality of the test is to make sure the fit works and
    we get a reasonable answer out **given** good initial guesses.
    """

    #
    # Create initial guesses - based on the s1_um_mJy_e1 
    # Multiplication by 1000 as the dispersion axis is 
    # automatically converted to AA
    #

    fg1 = models.Gaussian1D(amplitude=1900, mean=5575, stddev=0.02*1000,
            bounds = {'amplitude': (1700, 2100), 'mean': (5000, 6000), 'stddev': (1.1754943508222875e-38, 200)})

    fg2 = models.Gaussian1D(amplitude=400, mean=6232, stddev=0.01*1000,
            bounds = {'amplitude': (200, 600), 'mean': (5500, 7000), 'stddev': (1.1754943508222875e-38, 200)})

    fg3 = models.Gaussian1D(amplitude=-300, mean=8005, stddev=0.01*1000,
            bounds = {'amplitude': (-500, -100), 'mean': (7500, 8500), 'stddev': (1.1754943508222875e-38, 200)})

    framp = models.Linear1D(slope=0.8, intercept=0)

    #
    #  Fit the simulated spectrum based on the input models (including bounds)
    #

    fitted_models = fit_models(simulated_spectra.s1_um_mJy_e1, [fg1, fg2, fg3, framp])

    #
    # Given the output ``fitted_models`` (which fits the data to calculate amplitude
    # mean and stddev for each gaussian model and slope of the ramp), calculated
    # the resultant fitted flux.
    #

    fitted_flux = fitted_models(simulated_spectra.s1_um_mJy_e1.wavelength.value)

    #
    # These are not going to match as there is noise in the original data
    # but the RMS should come out to approximately 1285 (calculated manually
    # and after visually checking the fit).
    #

    rms = np.sqrt(np.sum((np.array(simulated_spectra.s1_um_mJy_e1.flux.value)-np.array(fitted_flux))**2))

    assert np.allclose(rms, 1285, atol=1)


def test_fitting_simple_levmarlsq(simulated_spectra):
    """
    This test fits the the first simulated spectrum from the fixture.  The
    initial guesses are manually set here with bounds that essentially make
    sense as the functionality of the test is to make sure the fit works and
    we get a reasonable answer out **given** good initial guesses.
    """

    #
    # Create initial guesses - based on the s1_um_mJy_e1 
    # Multiplication by 1000 as the dispersion axis is 
    # automatically converted to AA
    #

    fg1 = models.Gaussian1D(amplitude=1900, mean=5575, stddev=0.02*1000,
            bounds = {'amplitude': (1700, 2100), 'mean': (5000, 6000), 'stddev': (1.1754943508222875e-38, 200)})

    fg2 = models.Gaussian1D(amplitude=400, mean=6232, stddev=0.01*1000,
            bounds = {'amplitude': (200, 600), 'mean': (5500, 7000), 'stddev': (1.1754943508222875e-38, 200)})

    fg3 = models.Gaussian1D(amplitude=-300, mean=8005, stddev=0.01*1000,
            bounds = {'amplitude': (-500, -100), 'mean': (7500, 8500), 'stddev': (1.1754943508222875e-38, 200)})

    framp = models.Linear1D(slope=0.8, intercept=0)

    #
    #  Fit the simulated spectrum based on the input models (including bounds)
    #

    fitted_models = fit_models(simulated_spectra.s1_um_mJy_e1, [fg1, fg2, fg3, framp], 
                               fitter=fitting.LevMarLSQFitter)

    #
    # Given the output ``fitted_models`` (which fits the data to calculate amplitude
    # mean and stddev for each gaussian model and slope of the ramp), calculated
    # the resultant fitted flux.
    #

    fitted_flux = fitted_models(simulated_spectra.s1_um_mJy_e1.wavelength.value)

    #
    # These are not going to match as there is noise in the original data
    # but the RMS should come out to approximately 1285 (calculated manually
    # and after visually checking the fit).
    #

    rms = np.sqrt(np.sum((np.array(simulated_spectra.s1_um_mJy_e1.flux.value)-np.array(fitted_flux))**2))

    assert np.allclose(rms, 2715, atol=1)
