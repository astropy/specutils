import os
import pytest
import astropy.io.ascii as ascii
from astropy import units as u
import numpy as np


from specutils.io import read_fits


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_multispec_linear():
    iraf = ascii.read(data_path('multispec-linear-log.fits.0001.dat'),
                      Reader=ascii.NoHeader, names=['wave', 'flux'])
    spectra = read_fits.read_fits_spectrum1d(data_path('multispec-linear-log.fits'))
    spec = spectra[0]
    np.testing.assert_allclose(iraf['wave'], spec.dispersion.value)
    assert spec.dispersion.unit == u.Angstrom

def test_multispec_log_linear():
    iraf = ascii.read(data_path('log-linear.dat'), Reader=ascii.NoHeader, names = ['wave', 'flux'])
    spectra = read_fits.read_fits_spectrum1d(data_path('multispec-linear-log.fits'))
    spec = spectra[1]
    np.testing.assert_allclose(iraf['wave'], spec.dispersion.value)
    assert spec.dispersion.unit == u.Angstrom

def test_multispec_equispec_linear():
    iraf11 = ascii.read(data_path('multispec_equispec.11.dat'),
                      Reader=ascii.NoHeader, names=['wave', 'flux'])
    iraf12 = ascii.read(data_path('multispec_equispec.12.dat'),
                      Reader=ascii.NoHeader, names=['wave', 'flux'])
    iraf21 = ascii.read(data_path('multispec_equispec.21.dat'),
                      Reader=ascii.NoHeader, names=['wave', 'flux'])
    iraf22 = ascii.read(data_path('multispec_equispec.22.dat'),
                      Reader=ascii.NoHeader, names=['wave', 'flux'])
    spectra = read_fits.read_fits_spectrum1d(data_path('multispec_equispec.fits'))

    np.testing.assert_allclose(iraf11['wave'], spectra[0][0].dispersion.value)
    np.testing.assert_allclose(iraf12['wave'], spectra[0][1].dispersion.value)
    np.testing.assert_allclose(iraf21['wave'], spectra[1][0].dispersion.value)
    np.testing.assert_allclose(iraf22['wave'], spectra[1][1].dispersion.value)
    spec = spectra[0][0]
    assert spec.dispersion.unit == u.Angstrom

def test_1d_multispec_combined():
    legendre = read_fits.read_fits_spectrum1d(data_path('TRES.fits'))[0]
    combined, chebyshev = read_fits.read_fits_spectrum1d(
        data_path('Combined.fits'))[0:2]
    weight_legendre = 2.0
    offset_legendre = 1.0
    weight_chebyshev = 3.0
    offset_chebyshev = 0.0

    test_value = weight_legendre*(offset_legendre+legendre.dispersion.value)\
        + weight_chebyshev*(offset_chebyshev+chebyshev.dispersion.value)
    np.testing.assert_allclose(combined.dispersion.value, test_value)

def test_multispec_legendre():
    iraf = ascii.read(data_path('TRES.dat'), data_start = 127, Reader = ascii.NoHeader, names = ['wave', 'flux'])
    spectra = read_fits.read_fits_spectrum1d(data_path('TRES.fits'))
    spec = spectra[10]
    np.testing.assert_allclose(iraf['wave'], spec.dispersion.value)
    assert spec.dispersion.unit == u.Angstrom

def test_multispec_chebyshev():
    iraf = ascii.read(data_path('AAO_11.txt'), data_start = 175, Reader = ascii.NoHeader, names = ['wave', 'flux'])
    spectra = read_fits.read_fits_spectrum1d(data_path('AAO.fits'))
    spec = spectra[10]
    np.testing.assert_allclose(iraf['wave'], spec.wavelength.value)

def test_1dspec_UVES():
    spec = read_fits.read_fits_spectrum1d(data_path('UVES.fits'))
    iraf = ascii.read(data_path('uves_iraf_read_'
                                'truncated.dat'), names=['index', 'wave', 'flux'])
    np.testing.assert_allclose(spec.dispersion[iraf['index']], iraf['wave'])

    assert not hasattr(spec.dispersion, 'unit')

def test_1dspec_vrad():
    iraf = ascii.read(data_path('gbt_1d_iraf_read.dat'), names=['wave', 'flux'])
    spec = read_fits.read_fits_spectrum1d(data_path('gbt_1d.fits'))
    np.testing.assert_allclose(iraf['wave'], spec.dispersion.value)
    assert spec.dispersion.unit == u.Unit('km/s')

