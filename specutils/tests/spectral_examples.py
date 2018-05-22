import numpy as np
import astropy.units as u
from astropy.modeling import models
from specutils.spectra import Spectrum1D
import pytest

class SpectraExamples(object):

    def __init__(self):
        # Create the base wavelengths and flux
        self.wavelengths_um = np.linspace(0.4, 1.05, 100)

        g1 = models.Gaussian1D(amplitude=2000, mean=0.56, stddev=0.01)
        g2 = models.Gaussian1D(amplitude=500, mean=0.62, stddev=0.02)
        g3 = models.Gaussian1D(amplitude=-400, mean=0.80, stddev=0.02)
        g4 = models.Gaussian1D(amplitude=-350, mean=0.52, stddev=0.01)
        ramp = models.Linear1D(slope=300, intercept=0.0)

        self.base_flux = g1(self.wavelengths_um) + g2(self.wavelengths_um) + g3(self.wavelengths_um) +\
                         g4(self.wavelengths_um) + ramp(self.wavelengths_um) + 1000

        # Initialize the seed so the random numbers are not quite as random
        np.random.seed(42)

        # Create two spectra with the only difference in the instance of noise
        self._flux_e1 = self.base_flux + 400*np.random.random(self.base_flux.shape)
        self._s1_um_mJy_e1 = Spectrum1D(spectral_axis=self.wavelengths_um*u.um, flux=self._flux_e1*u.mJ)

        self._flux_e2 = self.base_flux + 400*np.random.random(self.base_flux.shape)
        self._s1_um_mJy_e2 = Spectrum1D(spectral_axis=self.wavelengths_um*u.um, flux=self._flux_e2*u.mJ)

        # Create on spectrum with the same flux but in angstrom units
        self.wavelengths_AA = self.wavelengths_um*10000
        self._s1_AA_mJy_e1 = Spectrum1D(spectral_axis=self.wavelengths_AA*u.AA, flux=self._flux_e1*u.mJ)

    # Create two examples of different noise
    @property
    def s1_um_mJy_e1_flux(self):
        return self._flux_e1

    @property
    def s1_um_mJy_e1(self):
        return self._s1_um_mJy_e1

    @property
    def s1_um_mJy_e2(self):
        return self._s1_um_mJy_e2

    @property
    def s1_um_mJy_e2_flux(self):
        return self._flux_e2

    @property
    def s1_AA_mJy_e1(self):
        return self._s1_AA_mJy_e1

    @property
    def s1_AA_mJy_e1_flux(self):
        return self._flux_e1


@pytest.fixture
def define_spectra():
    return SpectraExamples()
