import astropy.units as u
import numpy as np

from astropy.io import fits, misc
from astropy.utils.data import get_pkg_data_filename
from ..spectra.spectrum1d import Spectrum1D
from ..manipulation import align_spectra


def test_aligning():

    data = misc.fnunpickle(get_pkg_data_filename('data/mos-nrs1.pck'))

    wave1 = data['spectrum1']['wavelengths']
    spec1 = data['spectrum1']['flux']

    wave2 = data['spectrum2']['wavelengths']
    spec2 = data['spectrum2']['flux']

    spectrum1 = Spectrum1D(flux=spec1, spectral_axis=wave1)
    spectrum2 = Spectrum1D(flux=spec2, spectral_axis=wave2)

    out_spectrum = align_spectra(spectrum1, spectrum2)

    assert np.allclose(out_spectrum.spectral_axis.value, data['out_spectrum']['wavelengths'].value)
    assert out_spectrum.spectral_axis.unit == data['out_spectrum']['wavelengths'].unit

    assert np.allclose(out_spectrum.flux.value, data['out_spectrum']['flux'].value)
    assert out_spectrum.flux.unit == data['out_spectrum']['flux'].unit
