import numpy as np
import astropy.units as u
from astropy.modeling import models
from ..spectra import Spectrum1D
import pytest


class SpectraExamples(object):
    """
    The ``SpectralExamples`` class is a *container class* that has
    several examples of simple spectra that are to be used in the tests
    (e.g., arithmetic tests, smoothing tests etc).

    The purpose of this being a test class instead of using a `Spectrum1D`
    directly is that it contains both the `Spectrum1D` object and the flux
    that was used to *create* the Spectrum.  That's for tests that ensure
    the simpler operations just on the flux arrays are carried through to
    the `Spectrum1D` operations.

    Each of the spectra are created from a base noise-less spectrum
    constructed from 4 Gaussians and a ramp. Then three example spectra
    are created, and then gaussian random noise is added.

        1. s1_um_mJy_e1 - 4 Gaussians + ramp with one instantion of noise
                          dispersion: um, flux: mJy

        2. s1_um_mJy_e2 - same as 1, but with a different instance of noise
                          dispersion: um, flux: mJy

        3. s1_AA_mJy_e3 - same as 1, but with a third instance of noise
                          dispersion: Angstroms, flux: mJy

        4. s1_AA_nJy_e3 - same as 1, but with a fourth instance of noise
                          dispersion: Angstroms, flux: nJy
    """

    def __init__(self):

        #
        # Create the base wavelengths and flux
        #

        self.wavelengths_um = np.linspace(0.4, 1.05, 100)

        g1 = models.Gaussian1D(amplitude=2000, mean=0.56, stddev=0.01)
        g2 = models.Gaussian1D(amplitude=500, mean=0.62, stddev=0.02)
        g3 = models.Gaussian1D(amplitude=-400, mean=0.80, stddev=0.02)
        g4 = models.Gaussian1D(amplitude=-350, mean=0.52, stddev=0.01)
        ramp = models.Linear1D(slope=300, intercept=0.0)

        self.base_flux = g1(self.wavelengths_um) + g2(self.wavelengths_um) + \
                         g3(self.wavelengths_um) + g4(self.wavelengths_um) + \
                         ramp(self.wavelengths_um) + 1000

        #
        # Initialize the seed so the random numbers are not quite as random
        #

        np.random.seed(42)

        #
        # Create two spectra with the only difference in the instance of noise
        #

        self._flux_e1 = self.base_flux + 400 * np.random.random(self.base_flux.shape)
        self._s1_um_mJy_e1 = Spectrum1D(spectral_axis=self.wavelengths_um * u.um,
                                        flux=self._flux_e1 * u.mJy)

        self._flux_e2 = self.base_flux + 400 * np.random.random(self.base_flux.shape)
        self._s1_um_mJy_e2 = Spectrum1D(spectral_axis=self.wavelengths_um * u.um,
                                        flux=self._flux_e2 * u.mJy)

        #
        # Create on spectrum with the same flux but in angstrom units
        #

        self.wavelengths_AA = self.wavelengths_um * 10000
        self._s1_AA_mJy_e3 = Spectrum1D(spectral_axis=self.wavelengths_AA * u.AA,
                                        flux=self._flux_e1 * u.mJy)

        #
        # Create on spectrum with the same flux but in angstrom units and nJy
        #

        self._flux_e4 = (self.base_flux + 400 * np.random.random(self.base_flux.shape)) * 1000000
        self._s1_AA_nJy_e4 = Spectrum1D(spectral_axis=self.wavelengths_AA * u.AA,
                                        flux=self._flux_e4 * u.nJy)


    @property
    def s1_um_mJy_e1(self):
        return self._s1_um_mJy_e1

    @property
    def s1_um_mJy_e1_flux(self):
        return self._flux_e1

    @property
    def s1_um_mJy_e2(self):
        return self._s1_um_mJy_e2

    @property
    def s1_um_mJy_e2_flux(self):
        return self._flux_e2

    @property
    def s1_AA_mJy_e3(self):
        return self._s1_AA_mJy_e3

    @property
    def s1_AA_mJy_e3_flux(self):
        return self._flux_e1

    @property
    def s1_AA_nJy_e4(self):
        return self._s1_AA_nJy_e4

    @property
    def s1_AA_nJy_e4_flux(self):
        return self._flux_e4


@pytest.fixture
def simulated_spectra():
    """
    The method will be called as a fixture to tests.

    Parameters
    ----------
    N/A

    Return
    ------
    ``SpectralExamples``
        An instance of the SpectraExamples class.

    Examples
    --------
    This fixture can be used in a test as:

    ```
    from .spectral_examples import spectral_examples

    def test_add_spectra(spectral_examples):

        # Get the numpy array of data
        flux1 = define_spectra.s1_um_mJy_e1_flux
        flux2 = define_spectra.s1_um_mJy_e2_flux
        flux3 = flux1 + flux2

        # Calculate using the spectrum1d/nddata code
        spec3 = define_spectra.s1_um_mJy_e1 + define_spectra.s1_um_mJy_e2

        assert np.allclose(spec3.flux.value, flux3)
    ```

    """

    return SpectraExamples()
