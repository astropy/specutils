import astropy.units as u
import astropy.wcs as fitswcs
import gwcs
import numpy as np
import pytest
from astropy.nddata import StdDevUncertainty
from astropy.coordinates import SpectralCoord

from .conftest import remote_access
from ..spectra import Spectrum1D

def test_empty_spectrum():
    spec = Spectrum1D(spectral_axis=[]*u.um,
                      flux=[]*u.Jy)

    assert isinstance(spec.spectral_axis, SpectralCoord)
    assert spec.spectral_axis.size == 0

    assert isinstance(spec.flux, u.Quantity)
    assert spec.flux.size == 0


def test_create_from_arrays():
    spec = Spectrum1D(spectral_axis=np.arange(50) * u.AA,
                      flux=np.random.randn(50) * u.Jy)

    assert isinstance(spec.spectral_axis, SpectralCoord)
    assert spec.spectral_axis.size == 50

    assert isinstance(spec.flux, u.Quantity)
    assert spec.flux.size == 50

    # Test creating spectrum with unknown arguments
    with pytest.raises(ValueError) as e_info:
        spec = Spectrum1D(wavelength=np.arange(1, 50) * u.nm,
                          flux=np.random.randn(48) * u.Jy)


def test_create_from_multidimensional_arrays():
    """
    This is a test for a bug that was fixed by #283. It makes sure that
    multidimensional flux arrays are handled properly when creating Spectrum1D
    objects.
    """

    freqs = np.arange(50) * u.GHz
    flux = np.random.random((5, len(freqs))) * u.Jy
    spec = Spectrum1D(spectral_axis=freqs, flux=flux)

    assert (spec.frequency == freqs).all()
    assert (spec.flux == flux).all()

    # Mis-matched lengths should raise an exception (unless freqs is one longer
    # than flux, in which case it's interpreted as bin edges)
    freqs = np.arange(50) * u.GHz
    flux = np.random.random((5, len(freqs)-10)) * u.Jy
    with pytest.raises(ValueError) as e_info:
        spec = Spectrum1D(spectral_axis=freqs, flux=flux)

def test_create_from_quantities():
    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.randn(49) * u.Jy)

    assert isinstance(spec.spectral_axis, SpectralCoord)
    assert spec.spectral_axis.unit == u.nm
    assert spec.spectral_axis.size == 49

    # Mis-matched lengths should raise an exception (unless freqs is one longer
    # than flux, in which case it's interpreted as bin edges)
    with pytest.raises(ValueError) as e_info:
        spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.randn(47) * u.Jy)


def test_create_implicit_wcs():
    spec = Spectrum1D(spectral_axis=np.arange(50) * u.AA,
                      flux=np.random.randn(50) * u.Jy)

    assert isinstance(spec.wcs, gwcs.wcs.WCS)

    pix2world = spec.wcs.pixel_to_world(np.arange(5, 10))

    assert pix2world.size == 5
    assert isinstance(pix2world, np.ndarray)


def test_create_implicit_wcs_with_spectral_unit():
    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.randn(49) * u.Jy)

    assert isinstance(spec.wcs, gwcs.wcs.WCS)

    pix2world = spec.wcs.pixel_to_world(np.arange(5, 10))

    assert pix2world.size == 5
    assert isinstance(pix2world, np.ndarray)


def test_create_with_spectral_coord():

    spectral_coord = SpectralCoord(np.arange(5100, 5150)*u.AA, radial_velocity=u.Quantity(1000.0, "km/s"))
    flux = np.random.randn(50)*u.Jy
    spec = Spectrum1D(spectral_axis = spectral_coord, flux=flux)

    assert spec.radial_velocity == u.Quantity(1000.0, "km/s")
    assert isinstance(spec.spectral_axis, SpectralCoord)
    assert spec.spectral_axis.size == 50


def test_spectral_axis_conversions():
    # By default the spectral axis units should be set to angstroms
    spec = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy,
                      spectral_axis=np.array([400, 500]) * u.AA)

    assert np.all(spec.spectral_axis == np.array([400, 500]) * u.angstrom)
    assert spec.spectral_axis.unit == u.angstrom

    spec = Spectrum1D(spectral_axis=np.arange(50) * u.AA,
                      flux=np.random.randn(50) * u.Jy)

    assert spec.wavelength.unit == u.AA

    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.randn(49) * u.Jy)

    assert spec.frequency.unit == u.GHz

    with pytest.raises(ValueError) as e_info:
        spec.velocity

    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.randn(49) * u.Jy)

    new_spec = spec.with_spectral_unit(u.GHz)


@pytest.mark.parametrize('unit', ['micron', 'GHz', 'cm**-1', 'eV'])
def test_spectral_axis_equivalencies(unit):
    """Test that `u.spectral` equivalencies are enabled for `spectral_axis`."""

    spectral_axis=np.array([3400, 5000, 6660]) * u.AA
    spec = Spectrum1D(flux=np.array([26.0, 30.0, 44.5]) * u.Jy, spectral_axis=spectral_axis)

    new_axis = spectral_axis.to(unit, equivalencies=u.spectral())
    assert u.allclose(spec.spectral_axis.to(unit), new_axis)


def test_redshift():
    spec = Spectrum1D(flux=np.array([26.0, 30., 44.5]) * u.Jy,
                      spectral_axis=np.array([4000, 6000, 8000]) * u.AA,
                      velocity_convention='optical',
                      rest_value=6000 * u.AA)

    assert u.allclose(spec.velocity, [-99930.8, 0, 99930.8]*u.km/u.s,
                     atol=0.5*u.km/u.s)

    spec = Spectrum1D(flux=np.array([26.0, 30., 44.5]) * u.Jy,
                      spectral_axis=np.array([4000, 6000, 8000]) * u.AA,
                      velocity_convention='optical',
                      rest_value=6000 * u.AA,
                      redshift= 0.1)

    assert u.allclose(spec.velocity, [-71443.75318854, 28487.0661448, 128417.88547813]*u.km/u.s,
                      atol=0.5*u.km/u.s)

    #-------------------------

    spec = Spectrum1D(flux=np.array([26.0, 30.0, 44.5]) * u.Jy,
                      spectral_axis=np.array([10.5, 11.0, 11.5]) * u.GHz,
                      velocity_convention='radio',
                      rest_value=11.0 * u.GHz)

    assert u.allclose(spec.velocity, [13626., 0, -13626]*u.km/u.s,
                      atol=1*u.km/u.s)

    spec = Spectrum1D(flux=np.array([26.0, 30.0, 44.5]) * u.Jy,
                      spectral_axis=np.array([10.5, 11.0, 11.5]) * u.GHz,
                      velocity_convention='radio',
                      rest_value=11.0 * u.GHz,
                      redshift= 0.1)

    assert u.allclose(spec.velocity, [42113.99605389, 28487.0661448 , 14860.13623571]*u.km/u.s,
                      atol=1*u.km/u.s)

    #------------------------- radial velocity mode

    spec = Spectrum1D(flux=np.array([26.0, 30., 44.5]) * u.Jy,
                      spectral_axis=np.array([4000, 6000, 8000]) * u.AA,
                      velocity_convention='optical',
                      rest_value=6000 * u.AA)

    assert u.allclose(spec.velocity, [-99930.8, 0.0, 99930.8]*u.km/u.s,
                      atol=0.5*u.km/u.s)

    spec = Spectrum1D(flux=np.array([26.0, 30., 44.5]) * u.Jy,
                      spectral_axis=np.array([4000, 6000, 8000]) * u.AA,
                      velocity_convention='optical',
                      rest_value=6000 * u.AA,
                      radial_velocity=1000.*u.km/u.s)

    assert u.allclose(spec.velocity, [-98930.8, 1000.0, 100930.8]*u.km/u.s,
                      atol=0.5*u.km/u.s)


def test_flux_unit_conversion():
    # By default the flux units should be set to Jy
    s = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy,
                   spectral_axis=np.array([400, 500]) * u.nm)
    assert np.all(s.flux == np.array([26.0, 44.5]) * u.Jy)
    assert s.flux.unit == u.Jy

    # Simple Unit Conversion
    s = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy, spectral_axis=np.array([400, 500])*u.nm)
    converted_spec = s.new_flux_unit(unit=u.uJy)
    assert ((26.0 * u.Jy).to(u.uJy) == converted_spec.flux[0])

    # Make sure incompatible units raise UnitConversionError
    with pytest.raises(u.UnitConversionError):
        s.new_flux_unit(unit=u.m)

    # Pass custom equivalencies
    s = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy,
                   spectral_axis=np.array([400, 500]) * u.nm)
    eq = [[u.Jy, u.m,
          lambda x: np.full_like(np.array(x), 1000.0, dtype=np.double),
          lambda x: np.full_like(np.array(x), 0.001, dtype=np.double)]]
    converted_spec = s.new_flux_unit(unit=u.m, equivalencies=eq)
    assert 1000.0 * u.m == converted_spec.flux[0]

    # Check if suppressing the unit conversion works
    s = Spectrum1D(flux=np.array([26.0, 44.5]) * u.Jy, spectral_axis=np.array([400, 500]) * u.nm)
    new_spec = s.new_flux_unit("uJy", suppress_conversion=True)
    assert new_spec.flux[0] == 26.0 * u.uJy


def test_wcs_transformations():
    # Test with a GWCS
    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.randn(49) * u.Jy)

    pix_axis = spec.wcs.world_to_pixel(np.arange(20, 30) * u.nm)
    disp_axis = spec.wcs.pixel_to_world(np.arange(20, 30))

    assert isinstance(pix_axis, np.ndarray)
    assert isinstance(disp_axis, u.Quantity)

    # Test transform with different unit
    with u.set_enabled_equivalencies(u.spectral()):
        spec.wcs.world_to_pixel(np.arange(20, 30) * u.GHz)

    # Test with a FITS WCS
    my_wcs = fitswcs.WCS(header={'CDELT1': 1, 'CRVAL1': 6562.8, 'CUNIT1': 'Angstrom',
                                 'CTYPE1': 'WAVE', 'RESTFRQ': 1400000000, 'CRPIX1': 25},
                         naxis=1)

    spec = Spectrum1D(flux=[5,6,7] * u.Jy, wcs=my_wcs)

    pix_axis = spec.wcs.world_to_pixel(20 * u.um)
    disp_axis = spec.wcs.pixel_to_world(np.arange(20, 30))

    assert isinstance(pix_axis, np.ndarray)
    assert isinstance(disp_axis, u.Quantity)

    assert np.allclose(spec.wcs.world_to_pixel(7000*u.AA), [461.2])

def test_create_explicit_fitswcs():
    my_wcs = fitswcs.WCS(header={'CDELT1': 1, 'CRVAL1': 6562.8, 'CUNIT1': 'Angstrom',
                                 'CTYPE1': 'WAVE', 'RESTFRQ': 1400000000, 'CRPIX1': 25})

    spec = Spectrum1D(flux=[5,6,7] * u.Jy, wcs=my_wcs)
    spec = spec.with_velocity_convention("relativistic")

    assert isinstance(spec.spectral_axis, SpectralCoord)
    assert spec.spectral_axis.unit.is_equivalent(u.AA)

    pix2world = spec.wcs.pixel_to_world(np.arange(3))

    assert pix2world.size == 3

    assert isinstance(spec.wavelength, u.Quantity)
    assert spec.wavelength.size == 3
    assert spec.wavelength.unit == u.AA

    assert isinstance(spec.frequency, u.Quantity)
    assert spec.frequency.size == 3
    assert spec.frequency.unit == u.GHz

    assert isinstance(spec.velocity, u.Quantity)
    assert spec.velocity.size == 3
    assert spec.velocity.unit == u.Unit('km/s')


def test_create_with_uncertainty():
    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.sample(49) * u.Jy,
                      uncertainty=StdDevUncertainty(np.random.sample(49) * 0.1))

    assert isinstance(spec.uncertainty, StdDevUncertainty)

    spec = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm,
                      flux=np.random.sample(49) * u.Jy,
                      uncertainty=StdDevUncertainty(np.random.sample(49) * 0.1))

    assert spec.flux.unit == spec.uncertainty.unit

    # If flux and uncertainty are different sizes then raise exception
    wavelengths = np.arange(0, 10)
    flux=100*np.abs(np.random.randn(3, 4, 10))*u.Jy
    uncertainty = StdDevUncertainty(np.abs(np.random.randn(3, 2, 10))*u.Jy)

    with pytest.raises(ValueError) as e_info:
        s1d = Spectrum1D(spectral_axis=wavelengths*u.um,
                         flux=flux,
                         uncertainty=uncertainty)


@remote_access([{'id': '1481190', 'filename': 'L5g_0355+11_Cruz09.fits'}])
def test_read_linear_solution(remote_data_path):
    spec = Spectrum1D.read(remote_data_path, format='wcs1d-fits')

    assert isinstance(spec, Spectrum1D)

    assert isinstance(spec.flux, u.Quantity)
    assert isinstance(spec.spectral_axis, SpectralCoord)

    assert spec.flux.size == spec.data.size
    assert spec.spectral_axis.size == spec.data.size


def test_energy_photon_flux():
    spec = Spectrum1D(spectral_axis=np.linspace(100, 1000, 10) * u.nm,
                      flux=np.random.randn(10)*u.Jy)
    assert spec.energy.size == 10
    assert spec.photon_flux.size == 10
    assert spec.photon_flux.unit == u.photon * u.cm**-2 * u.s**-1 * u.nm**-1


def test_flux_nans_propagate_to_mask():
    """Check that indices in input flux with NaNs get propagated to the mask"""
    flux = np.random.randn(10)
    nan_idx = [0, 3, 5]
    flux[nan_idx] = np.nan
    spec = Spectrum1D(spectral_axis=np.linspace(100, 1000, 10) * u.nm,
                      flux=flux * u.Jy)
    assert spec.mask[nan_idx].all() == True


def test_repr():
    spec_with_wcs = Spectrum1D(spectral_axis=np.linspace(100, 1000, 10) * u.nm,
                               flux=np.random.random(10) * u.Jy)
    result = repr(spec_with_wcs)
    assert result.startswith('<Spectrum1D(flux=<Quantity [')
    assert 'spectral_axis=<SpectralAxis' in result

    spec_with_unc = Spectrum1D(spectral_axis=np.linspace(100, 1000, 10) * u.nm,
                               flux=np.random.random(10) * u.Jy,
                               uncertainty=StdDevUncertainty(
                                   np.random.sample(10), unit='Jy'))
    result = repr(spec_with_unc)
    assert result.startswith('<Spectrum1D(flux=<Quantity [')
    assert 'spectral_axis=<SpectralAxis' in result
    assert 'uncertainty=StdDevUncertainty(' in result


def test_str():
    spec = Spectrum1D(spectral_axis=np.linspace(100, 1000, 10) * u.nm,
                               flux=np.random.random(10) * u.Jy)
    result = str(spec)
    # Sanity check for contents of string representation
    assert result.startswith('Spectrum1D (length={})'.format(len(spec.flux)))
    lines = result.split('\n')
    flux = spec.flux
    sa = spec.spectral_axis
    assert len(lines) == 3
    assert lines[1].startswith('flux:')
    assert '[ {:.5}, ..., {:.5} ]'.format(flux[0], flux[-1]) in lines[1]
    assert 'mean={:.5}'.format(np.mean(flux)) in lines[1]
    assert lines[2].startswith('spectral axis:')
    assert '[ {:.5}, ..., {:.5} ]'.format(sa[0], sa[-1]) in lines[2]
    assert 'mean={:5}'.format(np.mean(sa)) in lines[2]

    # Test string representation with uncertainty
    spec_with_uncertainty = Spectrum1D(spectral_axis=np.linspace(100, 1000, 10) * u.nm,
                                       flux=np.random.random(10) * u.Jy,
                                       uncertainty=np.random.random(10))
    result = str(spec_with_uncertainty)
    lines = result.split('\n')
    unc = spec_with_uncertainty.uncertainty
    print(spec_with_uncertainty)
    assert len(lines) == 4
    assert lines[3].startswith('uncertainty')
    assert '[ {}, ..., {} ]'.format(unc[0], unc[-1]) in lines[3]

    # Test string representation with multiple flux
    spec_multi_flux = Spectrum1D(spectral_axis=np.linspace(100, 1000, 10) * u.nm,
                                 flux=np.random.random((3,10)) * u.Jy)
    result = str(spec_multi_flux)
    lines = result.split('\n')
    assert len(lines) == 5
    for i, line in enumerate(lines[1:4]):
        assert line.startswith('flux{:2}:'.format(i))

    # Test string representation with single-dimensional flux
    spec_single_flux = Spectrum1D([1] * u.Jy, [0] * u.nm)
    result = str(spec_single_flux)
    assert result == \
"""Spectrum1D (length=1)
flux:             [ 1.0 Jy ],  mean=1.0 Jy
spectral axis:    [ 0.0 nm ],  mean=0.0 nm"""


def test_equivalencies():
    """
    Test that after import `u.spectral` equivalencies are not enabled in the global namespace.
    """
    assert u.micron.is_equivalent(u.cm**-1) is False
    assert u.micron.is_equivalent(u.Hz) is False
    assert u.micron.is_equivalent(u.eV) is False
    assert u.Hz.is_equivalent(u.cm**-1) is False
    assert u.Hz.is_equivalent(u.eV) is False
