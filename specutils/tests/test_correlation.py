import astropy.units as u
import numpy as np
from astropy.nddata import StdDevUncertainty

from ..spectra.spectrum1d import Spectrum1D
from ..analysis import correlation
from astropy.tests.helper import quantity_allclose

SIZE = 41


def test_correlation():
    """
    Test correlation when both observed and template spectra have the same wavelength axis
    """
    # Seed np.random so that results are consistent
    np.random.seed(41)

    # Create test spectra
    spec_axis = np.linspace(0, SIZE, SIZE) * u.AA
    spec = Spectrum1D(spectral_axis=spec_axis,
                      flux=np.random.randn(SIZE) * u.Jy,
                      uncertainty=StdDevUncertainty(np.random.sample(SIZE), unit='Jy'))

    spec1 = Spectrum1D(spectral_axis=spec_axis,
                       flux=np.random.randn(SIZE) * u.Jy,
                       uncertainty=StdDevUncertainty(np.random.sample(SIZE), unit='Jy'))

    # Get result from correlation
    corr_result = correlation.template_correlate(spec, spec1)

    # Checks that lag at mid-point is zero.
    assert int((corr_result.spectral_axis[int(len(corr_result.spectral_axis)/2)]).value) == 0


