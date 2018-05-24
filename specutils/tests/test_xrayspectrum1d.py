import numpy as np
from scipy import optimize
import astropy.units as u
from astropy.modeling.powerlaws import PowerLaw1D
from specutils import XraySpectrum1D

def test_create_from_arrays():
    # Test that XraySpectrum1D can be initialized
    amp0, alpha0 = 3.e-3, 2.0
    powlaw0 = PowerLaw1D(amplitude=amp0, alpha=alpha0, x_0=1.e3)
    energy = np.linspace(0.2, 10.0, 8000)
    elo = energy[:-1]
    ehi = energy[1:]
    emid = 0.5 * (elo + ehi)
    counts = np.random.poisson(lam=powlaw0(emid), size=len(emid))
    test_spec = XraySpectrum1D(elo, ehi, u.keV, counts, exposure=1.0)
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
    new_model = PowerLaw1D(amplitude=5.e-3, alpha=1.8, x_0=1.e3)
    model_flux = new_model(test_spec.spectral_axis.value)
    ymodel = test_spec.apply_resp(model_flux)
    assert len(ymodel) == len(test_spec.spectral_axis)
    return
