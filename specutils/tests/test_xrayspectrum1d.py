import numpy as np
import astropy.units as u
from astropy.modeling.powerlaws import PowerLaw1D
from specutils import XraySpectrum1D, ARF, RMF

def test_create_from_arrays():
    # Test that XraySpectrum1D can be initialized
    amp0, alpha0 = 3.e-3, 2.0
    powlaw0 = PowerLaw1D(amplitude=amp0, alpha=alpha0, x_0=1.e3)
    energy = np.linspace(0.2, 10.0, 8000)
    elo = energy[:-1] * u.keV
    ehi = energy[1:] * u.keV
    emid = 0.5 * (elo + ehi)
    counts = np.random.poisson(lam=powlaw0(emid.value), size=len(emid)) * u.ct
    test_spec = XraySpectrum1D(elo, ehi, counts, exposure=1.0*u.second)
    return test_spec

def test_call_spectrum():
    # Test that the XraySpectrum1D has attributes that make it unique
    # to X-ray spectrum
    test_spec = test_create_from_arrays()
    assert len(test_spec.counts) == len(test_spec.spectral_axis)
    assert len(test_spec.bin_lo) == len(test_spec.spectral_axis)
    assert len(test_spec.bin_hi) == len(test_spec.spectral_axis)
    return

def test_apply_model():
    # Test that one can evaluate XraySpectrum1D with a model
    test_spec = test_create_from_arrays()
    new_model = PowerLaw1D(amplitude=5.e-3, alpha=1.8, x_0=1.0)
    model_flux = new_model(test_spec.spectral_axis.value)
    ymodel = test_spec.apply_response(model_flux)
    assert len(ymodel) == len(test_spec.spectral_axis)
    # Since there's no arf or rmf, should get the same array out
    assert np.all(ymodel == model_flux)

    # Now add an arf
    specresp = np.ones_like(test_spec.bin_lo.value) * 2.0
    arf = ARF(test_spec.bin_lo, test_spec.bin_hi,
                       specresp, test_spec.exposure)
    test_spec.assign_arf(arf)

    # Test that the arf is applied properly
    ymodel2 = test_spec.apply_response(model_flux)
    # The arf also gets multiplied by the exposure time, so that needs
    # to be included in this check
    check_value = ymodel2 / (model_flux * test_spec.exposure)
    assert np.all(check_value == 2.0)
    return

## Currently there is no RMF check because it requires a test data file
## However, we can test that it can be initialized with empty values
def test_init_rmf():
    rmf = RMF()
    return
