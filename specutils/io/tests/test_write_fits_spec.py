import os
from specutils.io import read_fits, write_fits
import numpy as np
import pytest

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_linear_write1():
    spec = read_fits.read_fits_spectrum1d(data_path('UVES.fits'))
    write_fits.write(spec, 'test.fits')
    test_spec = read_fits.read_fits_spectrum1d('test.fits')
    np.testing.assert_allclose(spec.data, test_spec.data)
    np.testing.assert_allclose(spec.dispersion, test_spec.dispersion)
    if hasattr(spec.dispersion, "unit") or hasattr(test_spec.dispersion, "unit"):
        assert spec.dispersion.unit == test_spec.dispersion.unit

def test_linear_write2():
    spec = read_fits.read_fits_spectrum1d(data_path('gbt_1d.fits'))
    write_fits.write(spec, 'test.fits')
    test_spec = read_fits.read_fits_spectrum1d('test.fits')
    np.testing.assert_allclose(spec.data, test_spec.data)
    np.testing.assert_allclose(spec.dispersion, test_spec.dispersion)
    if hasattr(spec.dispersion, "unit") or hasattr(test_spec.dispersion, "unit"):
        assert spec.dispersion.unit == test_spec.dispersion.unit

def test_multispec_legendre():
    spectra = read_fits.read_fits_spectrum1d(data_path('TRES.fits'))
    write_fits.write(spectra, 'test.fits')
    test_spectra = read_fits.read_fits_spectrum1d('test.fits')
    for spec, test_spec in zip(spectra, test_spectra):
        np.testing.assert_allclose(spec.data, test_spec.data)
        np.testing.assert_allclose(spec.dispersion, test_spec.dispersion)
    if hasattr(spec.dispersion, "unit") or hasattr(test_spec.dispersion, "unit"):
        assert spec.dispersion.unit == test_spec.dispersion.unit

def test_multispec_chebyshev():
    spectra = read_fits.read_fits_spectrum1d(data_path('AAO.fits'))
    write_fits.write(spectra, 'test.fits')
    test_spectra = read_fits.read_fits_spectrum1d('test.fits')
    for spec, test_spec in zip(spectra, test_spectra):
        np.testing.assert_allclose(spec.data, test_spec.data)
        np.testing.assert_allclose(spec.dispersion, test_spec.dispersion)
    if hasattr(spec.dispersion, "unit") or hasattr(test_spec.dispersion, "unit"):
        assert spec.dispersion.unit == test_spec.dispersion.unit

def test_multispec_combined():
    spectra = read_fits.read_fits_spectrum1d(data_path('Combined.fits'))
    write_fits.write(spectra, 'test.fits')
    test_spectra = read_fits.read_fits_spectrum1d('test.fits')
    for spec, test_spec in zip(spectra, test_spectra):
        np.testing.assert_allclose(spec.data, test_spec.data)
        np.testing.assert_allclose(spec.dispersion, test_spec.dispersion)
    if hasattr(spec.dispersion, "unit") or hasattr(test_spec.dispersion, "unit"):
        assert spec.dispersion.unit == test_spec.dispersion.unit

def test_multispec_spline():
    pytest.importorskip("scipy")
    spectra = read_fits.read_fits_spectrum1d(data_path('spline.fits'))
    write_fits.write(spectra, 'test.fits')
    test_spectra = read_fits.read_fits_spectrum1d('test.fits')
    for spec, test_spec in zip(spectra, test_spectra):
        np.testing.assert_allclose(spec.data, test_spec.data)
        np.testing.assert_allclose(spec.dispersion, test_spec.dispersion)
    if hasattr(spec.dispersion, "unit") or hasattr(test_spec.dispersion, "unit"):
        assert spec.dispersion.unit == test_spec.dispersion.unit

def test_multispec_linear():
    spectra = read_fits.read_fits_spectrum1d(data_path(
                                                'multispec-linear-log.fits'))
    write_fits.write(spectra, 'test.fits')
    test_spectra = read_fits.read_fits_spectrum1d('test.fits')
    for spec, test_spec in zip(spectra, test_spectra):
        np.testing.assert_allclose(spec.data, test_spec.data)
        np.testing.assert_allclose(spec.dispersion, test_spec.dispersion)
    if hasattr(spec.dispersion, "unit") or hasattr(test_spec.dispersion, "unit"):
        assert spec.dispersion.unit == test_spec.dispersion.unit

# TODO: create tests that utilize all methods to write fits header