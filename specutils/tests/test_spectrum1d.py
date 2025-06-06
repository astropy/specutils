import astropy.units as u
import astropy.wcs as fitswcs
import gwcs
import numpy as np
from numpy.testing import assert_allclose
import pytest

from astropy.nddata import StdDevUncertainty
from astropy.coordinates import SpectralCoord
from astropy.tests.helper import quantity_allclose
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.wcs import WCS

from .conftest import remote_access
from ..spectra import Spectrum1D, Spectrum


def test_empty_spectrum():
    spec = Spectrum(spectral_axis=[]*u.um,
                      flux=[]*u.Jy)

    assert isinstance(spec.spectral_axis, SpectralCoord)
    assert spec.spectral_axis.size == 0

    assert isinstance(spec.flux, u.Quantity)
    assert spec.flux.size == 0


def test_create_from_arrays():
    spec = Spectrum(spectral_axis=np.arange(50) * u.AA,
                    flux=np.ones(50) * u.Jy)

    assert isinstance(spec.spectral_axis, SpectralCoord)
    assert spec.spectral_axis.size == 50

    assert isinstance(spec.flux, u.Quantity)
    assert spec.flux.size == 50

    # Test creating spectrum with unknown arguments
    with pytest.raises(ValueError):
        spec = Spectrum(wavelength=np.arange(1, 50) * u.nm,
                        flux=np.ones(48) * u.Jy)


def test_create_from_multidimensional_arrays():
    """
    This is a test for a bug that was fixed by #283. It makes sure that
    multidimensional flux arrays are handled properly when creating Spectrum
    objects.
    """

    freqs = np.arange(50) * u.GHz
    flux = np.ones((5, len(freqs))) * u.Jy
    spec = Spectrum(spectral_axis=freqs, flux=flux)

    assert (spec.frequency == freqs).all()
    assert (spec.flux == flux).all()

    # Mis-matched lengths should raise an exception (unless freqs is one longer
    # than flux, in which case it's interpreted as bin edges)
    flux = np.ones((5, len(freqs) - 10)) * u.Jy
    with pytest.raises(ValueError):
        spec = Spectrum(spectral_axis=freqs, flux=flux)


def test_create_from_quantities():
    wav = np.arange(1, 50) * u.nm
    flux = np.ones(49) * u.Jy
    spec = Spectrum(spectral_axis=wav, flux=flux)

    assert isinstance(spec.spectral_axis, SpectralCoord)
    assert spec.spectral_axis.unit == u.nm
    assert spec.spectral_axis.size == 49

    # Mis-matched lengths should raise an exception (unless freqs is one longer
    # than flux, in which case it's interpreted as bin edges)
    with pytest.raises(ValueError):
        spec = Spectrum(spectral_axis=wav, flux=np.ones(47) * u.Jy)


def test_create_implicit_wcs():
    spec = Spectrum(spectral_axis=np.arange(50) * u.AA,
                    flux=np.ones(50) * u.Jy)

    assert isinstance(spec.wcs, gwcs.wcs.WCS)

    pix2world = spec.wcs.pixel_to_world(np.arange(5, 10))

    assert pix2world.size == 5
    assert isinstance(pix2world, np.ndarray)


def test_create_implicit_wcs_with_spectral_unit():
    spec = Spectrum(spectral_axis=np.arange(1, 50) * u.nm,
                    flux=np.ones(49) * u.Jy)

    assert isinstance(spec.wcs, gwcs.wcs.WCS)

    pix2world = spec.wcs.pixel_to_world(np.arange(5, 10))

    assert pix2world.size == 5
    assert isinstance(pix2world, np.ndarray)


def test_create_with_spectral_coord():

    spectral_coord = SpectralCoord(np.arange(5100, 5150) * u.AA, radial_velocity=u.Quantity(1000.0, "km/s"))
    flux = np.ones(50) * u.Jy
    spec = Spectrum(spectral_axis=spectral_coord, flux=flux)

    assert spec.radial_velocity == u.Quantity(1000.0, "km/s")
    assert isinstance(spec.spectral_axis, SpectralCoord)
    assert spec.spectral_axis.size == 50


def test_create_from_cube():

    flux = np.arange(24).reshape([2,3,4])*u.Jy
    wcs_dict = {"CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN", "CTYPE3": "WAVE-LOG",
                "CRVAL1": 205, "CRVAL2": 27, "CRVAL3": 3.622e-7,
                "CDELT1": -0.0001, "CDELT2": 0.0001, "CDELT3": 8e-11,
                "CRPIX1": 0, "CRPIX2": 0, "CRPIX3": 0}
    w = WCS(wcs_dict)

    spec = Spectrum(flux=flux, wcs=w)
    spec_axis_from_wcs = (np.exp(np.array([1,2])*w.wcs.cdelt[-1]/w.wcs.crval[-1]) *
                          w.wcs.crval[-1]*spec.spectral_axis.unit)

    assert spec.flux.shape == (2,3,4)
    assert spec.flux[1,2,3] == 23*u.Jy
    assert quantity_allclose(spec.spectral_axis, spec_axis_from_wcs)

    with pytest.raises(ValueError):
        spec2 = Spectrum(flux=flux, wcs=w, move_spectral_axis='Bad string')

    # Test moving spectral axis from first to last
    spec2 = Spectrum(flux=flux, wcs=w, move_spectral_axis='last')
    assert spec2.flux.shape == (4,3,2)
    assert spec2.flux[3,2,1] == 23*u.Jy
    assert quantity_allclose(spec2.spectral_axis, spec_axis_from_wcs)

    # Test moving spectral axis from last to first
    spec3 = Spectrum(flux=spec2.flux, wcs=spec2.wcs, move_spectral_axis='first')
    assert spec3.flux.shape == (2,3,4)
    assert spec3.flux[1,2,3] == 23*u.Jy
    assert quantity_allclose(spec3.spectral_axis, spec_axis_from_wcs)


def test_create_old_classname():
    """
    This is a test for the Spectrum subclass introduced in 1.20 to ease transition
    to the new class name in 2.0.
    """

    freqs = np.arange(50) * u.GHz
    flux = np.ones((5, len(freqs))) * u.Jy
    # The move_spectral_axis and spectral_axis_index keywords sre simply ignored for now.
    with pytest.warns(AstropyDeprecationWarning):
        spec = Spectrum1D(spectral_axis=freqs, flux=flux, spectral_axis_index=-1)
    assert_allclose(spec.flux, flux)
    assert_allclose(spec.spectral_axis.value, freqs.value)


def test_spectral_axis_conversions():
    # By default the spectral axis units should be set to angstroms
    spec = Spectrum(flux=np.array([26.0, 44.5]) * u.Jy,
                      spectral_axis=np.array([400, 500]) * u.AA)

    assert np.all(spec.spectral_axis == np.array([400, 500]) * u.angstrom)
    assert spec.spectral_axis.unit == u.angstrom

    flux = np.ones(49) * u.Jy
    spec = Spectrum(spectral_axis=np.arange(flux.size) * u.AA, flux=flux)

    assert spec.wavelength.unit == u.AA

    spec = Spectrum(spectral_axis=np.arange(1, 50) * u.nm, flux=flux)

    assert spec.frequency.unit == u.GHz

    with pytest.raises(ValueError):
        spec.velocity

    spec = Spectrum(spectral_axis=np.arange(100, 150) * u.nm, flux=flux)

    new_spec = spec.with_spectral_axis_unit(u.km / u.s, rest_value=125 * u.um,
                                            velocity_convention="relativistic")

    assert new_spec.spectral_axis.unit == u.km / u.s
    assert new_spec.wcs.world_axis_units[0] == "km.s**-1"
    # Make sure meta stored the old WCS correctly
    assert new_spec.meta["original_wcs"].world_axis_units[0] == "nm"
    assert new_spec.meta["original_spectral_axis_unit"] == "nm"

    wcs_dict = {"CTYPE1": "WAVE", "CRVAL1": 3.622e3, "CDELT1": 8e-2,
                "CRPIX1": 0, "CUNIT1": "Angstrom"}
    wcs_spec = Spectrum(flux=flux, wcs=WCS(wcs_dict),
                          meta={'header': wcs_dict.copy()})
    new_spec = wcs_spec.with_spectral_axis_unit(u.km / u.s, rest_value=125 * u.um,
                                                velocity_convention="relativistic")
    new_spec.meta['original_wcs'].wcs.crval = [3.777e-7]
    new_spec.meta['header']['CRVAL1'] = 3777.0

    assert wcs_spec.wcs.wcs.crval[0] == 3.622e-7
    assert wcs_spec.meta['header']['CRVAL1'] == 3622.


def test_spectral_axis_and_flux_conversions():
    """A little bit from both sets of tests."""
    spec = Spectrum(spectral_axis=np.arange(100, 150) * u.nm,
                      flux=np.ones(49) * u.Jy)

    new_spec = spec.with_spectral_axis_and_flux_units(
        u.km / u.s, u.uJy, rest_value=125 * u.um, velocity_convention="relativistic")

    assert new_spec.spectral_axis.unit == u.km/u.s
    assert new_spec.wcs.world_axis_units[0] == "km.s**-1"
    # Make sure meta stored the old WCS correctly
    assert new_spec.meta["original_wcs"].world_axis_units[0] == "nm"
    assert new_spec.meta["original_spectral_axis_unit"] == "nm"
    assert new_spec.flux.unit == u.uJy
    assert_allclose(new_spec.flux.value, 1000000)


def test_spectral_slice():
    spec = Spectrum(spectral_axis=np.linspace(100, 1000, 10) * u.nm,
                    flux=np.ones(10) * u.Jy)
    sliced_spec = spec[300*u.nm:600*u.nm]
    assert np.all(sliced_spec.spectral_axis == [300, 400, 500] * u.nm)

    sliced_spec = spec[300*u.nm:605*u.nm]
    assert np.all(sliced_spec.spectral_axis == [300, 400, 500, 600] * u.nm)

    sliced_spec = spec[:300*u.nm]
    assert np.all(sliced_spec.spectral_axis == [100, 200] * u.nm)

    sliced_spec = spec[800*u.nm:]
    assert np.all(sliced_spec.spectral_axis == [800, 900, 1000] * u.nm)

    # Test higher dimensional slicing
    spec = Spectrum(spectral_axis=np.linspace(100, 1000, 10) * u.nm,
                    flux=np.ones((10, 10)) * u.Jy,
                       spectral_axis_index=1)
    sliced_spec = spec[300*u.nm:600*u.nm]
    assert np.all(sliced_spec.spectral_axis == [300, 400, 500] * u.nm)

    sliced_spec = spec[4:6, 300*u.nm:600*u.nm]
    assert sliced_spec.shape == (2, 3)


@pytest.mark.parametrize('unit', ['micron', 'GHz', 'cm**-1', 'eV'])
def test_spectral_axis_equivalencies(unit):
    """Test that `u.spectral` equivalencies are enabled for `spectral_axis`."""

    spectral_axis=np.array([3400, 5000, 6660]) * u.AA
    spec = Spectrum(flux=np.array([26.0, 30.0, 44.5]) * u.Jy, spectral_axis=spectral_axis)

    new_axis = spectral_axis.to(unit, equivalencies=u.spectral())
    assert u.allclose(spec.spectral_axis.to(unit), new_axis)


def test_redshift():
    spec = Spectrum(flux=np.array([26.0, 30., 44.5]) * u.Jy,
                      spectral_axis=np.array([4000, 6000, 8000]) * u.AA,
                      velocity_convention='optical',
                      rest_value=6000 * u.AA)

    assert u.allclose(spec.velocity, [-99930.8, 0, 99930.8]*u.km/u.s,
                     atol=0.5*u.km/u.s)

    spec = Spectrum(flux=np.array([26.0, 30., 44.5]) * u.Jy,
                      spectral_axis=np.array([4000, 6000, 8000]) * u.AA,
                      velocity_convention='optical',
                      rest_value=6000 * u.AA,
                      redshift= 0.1)

    assert u.allclose(spec.velocity, [-71443.75318854, 28487.0661448, 128417.88547813]*u.km/u.s,
                      atol=0.5*u.km/u.s)

    # -------------------------

    spec = Spectrum(flux=np.array([26.0, 30.0, 44.5]) * u.Jy,
                      spectral_axis=np.array([10.5, 11.0, 11.5]) * u.GHz,
                      velocity_convention='radio',
                      rest_value=11.0 * u.GHz)

    assert u.allclose(spec.velocity, [13626., 0, -13626]*u.km/u.s,
                      atol=1*u.km/u.s)

    spec = Spectrum(flux=np.array([26.0, 30.0, 44.5]) * u.Jy,
                      spectral_axis=np.array([10.5, 11.0, 11.5]) * u.GHz,
                      velocity_convention='radio',
                      rest_value=11.0 * u.GHz,
                      redshift= 0.1)

    assert u.allclose(spec.velocity, [42113.99605389, 28487.0661448 , 14860.13623571]*u.km/u.s,
                      atol=1*u.km/u.s)

    # ------------------------- radial velocity mode

    spec = Spectrum(flux=np.array([26.0, 30., 44.5]) * u.Jy,
                      spectral_axis=np.array([4000, 6000, 8000]) * u.AA,
                      velocity_convention='optical',
                      rest_value=6000 * u.AA)

    assert u.allclose(spec.velocity, [-99930.8, 0.0, 99930.8]*u.km/u.s,
                      atol=0.5*u.km/u.s)

    spec = Spectrum(flux=np.array([26.0, 30., 44.5]) * u.Jy,
                      spectral_axis=np.array([4000, 6000, 8000]) * u.AA,
                      velocity_convention='optical',
                      rest_value=6000 * u.AA,
                      radial_velocity=1000.*u.km/u.s)

    assert u.allclose(spec.velocity, [-98930.8, 1000.0, 100930.8]*u.km/u.s,
                      atol=0.5*u.km/u.s)


def test_flux_unit_conversion():
    # By default the flux units should be set to Jy
    s = Spectrum(flux=np.array([26.0, 44.5]) * u.Jy,
                   spectral_axis=np.array([400, 500]) * u.nm)
    assert np.all(s.flux == np.array([26.0, 44.5]) * u.Jy)
    assert s.flux.unit == u.Jy

    # Simple Unit Conversion
    s = Spectrum(flux=np.array([26.0, 44.5]) * u.Jy, spectral_axis=np.array([400, 500])*u.nm)
    converted_spec = s.with_flux_unit(unit=u.uJy)
    assert ((26.0 * u.Jy).to(u.uJy) == converted_spec.flux[0])

    # Make sure incompatible units raise UnitConversionError
    with pytest.raises(u.UnitConversionError):
        s.with_flux_unit(unit=u.m)

    # Pass custom equivalencies
    s = Spectrum(flux=np.array([26.0, 44.5]) * u.Jy,
                   spectral_axis=np.array([400, 500]) * u.nm)
    eq = [[u.Jy, u.m,
          lambda x: np.full_like(np.array(x), 1000.0, dtype=np.double),
          lambda x: np.full_like(np.array(x), 0.001, dtype=np.double)]]
    converted_spec = s.with_flux_unit(unit=u.m, equivalencies=eq)
    assert 1000.0 * u.m == converted_spec.flux[0]

    # Check if suppressing the unit conversion works
    s = Spectrum(flux=np.array([26.0, 44.5]) * u.Jy, spectral_axis=np.array([400, 500]) * u.nm)
    new_spec = s.with_flux_unit("uJy", suppress_conversion=True)
    assert new_spec.flux[0] == 26.0 * u.uJy


def test_wcs_transformations():
    # Test with a GWCS
    spec = Spectrum(spectral_axis=np.arange(1, 50) * u.nm,
                    flux=np.ones(49) * u.Jy)

    # After spacetelescope/gwcs#457 is merged and released, this can be changed to
    # pix_axis = spec.wcs.world_to_pixel_values(np.arange(20, 30) * u.nm)
    pix_axis = spec.wcs.world_to_pixel(SpectralCoord(np.arange(20, 30) * u.nm))
    disp_axis = spec.wcs.pixel_to_world(np.arange(20, 30))

    assert isinstance(pix_axis, np.ndarray)
    assert isinstance(disp_axis, u.Quantity)

    # Test transform with different unit
    with u.set_enabled_equivalencies(u.spectral()):
        # After spacetelescope/gwcs#457 is merged and released, this can be changed to
        # spec.wcs.world_to_pixel_values(np.arange(20, 30) * u.GHz)
        spec.wcs.world_to_pixel(SpectralCoord(np.arange(6.9e6, 1e7) * u.GHz))

    # Test with a FITS WCS
    my_wcs = fitswcs.WCS(header={'CDELT1': 1, 'CRVAL1': 6562.8, 'CUNIT1': 'Angstrom',
                                 'CTYPE1': 'WAVE', 'RESTFRQ': 1400000000, 'CRPIX1': 25},
                         naxis=1)

    spec = Spectrum(flux=[5,6,7] * u.Jy, wcs=my_wcs)

    pix_axis = spec.wcs.world_to_pixel(20 * u.um)
    disp_axis = spec.wcs.pixel_to_world(np.arange(20, 30))

    assert isinstance(pix_axis, np.ndarray)
    assert isinstance(disp_axis, u.Quantity)

    assert np.allclose(spec.wcs.world_to_pixel(7000*u.AA), [461.2])


def test_create_explicit_fitswcs():
    my_wcs = fitswcs.WCS(header={'CDELT1': 1, 'CRVAL1': 6562.8, 'CUNIT1': 'Angstrom',
                                 'CTYPE1': 'WAVE', 'RESTFRQ': 1400000000, 'CRPIX1': 25})

    spec = Spectrum(flux=[5,6,7] * u.Jy, wcs=my_wcs)
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
    spec = Spectrum(spectral_axis=np.arange(1, 50) * u.nm,
                    flux=np.ones(49) * u.Jy,
                    uncertainty=StdDevUncertainty(np.ones(49) * 0.1))

    assert isinstance(spec.uncertainty, StdDevUncertainty)
    assert spec.flux.unit == spec.uncertainty.unit

    # If flux and uncertainty are different sizes then raise exception
    wavelengths = np.arange(10) * u.um
    flux= np.ones((3, 4, 10)) * u.Jy
    uncertainty = StdDevUncertainty(np.ones((3, 2, 10)) * u.Jy)

    with pytest.raises(ValueError):
        Spectrum(spectral_axis=wavelengths, flux=flux, uncertainty=uncertainty)


@pytest.mark.parametrize("flux_unit", ["adu", "ct/s", "count"])
def test_flux_unit_io_roundtrip(tmp_path, flux_unit):
    # regression test for https://github.com/astropy/specutils/pull/1018
    fname = str(tmp_path / 'flux_unit_io_roundtrip.fits')
    sp = Spectrum(flux=np.ones(11) * u.Unit(flux_unit),
                    spectral_axis=np.arange(1, 12) * u.Unit('Hz'))
    sp.write(fname, overwrite=True)

    sp_load = Spectrum.read(fname)
    assert sp_load.flux.unit == sp.flux.unit


@pytest.mark.filterwarnings('ignore::astropy.io.fits.verify.VerifyWarning')
@remote_access([{'id': '1481190', 'filename': 'L5g_0355+11_Cruz09.fits'}])
def test_read_linear_solution(remote_data_path):
    spec = Spectrum.read(remote_data_path, format='wcs1d-fits')

    assert isinstance(spec, Spectrum)

    assert isinstance(spec.flux, u.Quantity)
    assert isinstance(spec.spectral_axis, SpectralCoord)

    assert spec.flux.size == spec.data.size
    assert spec.spectral_axis.size == spec.data.size


def test_energy_photon_flux():
    spec = Spectrum(spectral_axis=np.linspace(100, 1000, 10) * u.nm,
                    flux=np.ones(10) * u.Jy)
    assert spec.energy.size == 10
    assert spec.photon_flux.size == 10
    assert spec.photon_flux.unit == u.photon * u.cm**-2 * u.s**-1 * u.nm**-1


def test_flux_nans_propagate_to_mask():
    """Check that indices in input flux with NaNs get propagated to the mask"""
    flux = np.ones(10)
    nan_idx = [0, 3, 5]
    flux[nan_idx] = np.nan
    spec = Spectrum(spectral_axis=np.linspace(100, 1000, 10) * u.nm,
                      flux=flux * u.Jy)
    assert spec.mask[nan_idx].all() == True  # noqa


def test_repr():
    wav = np.linspace(100, 1000, 10) * u.nm
    flux = np.ones(10) * u.Jy
    spec_with_wcs = Spectrum(spectral_axis=wav, flux=flux)
    result = repr(spec_with_wcs)
    assert result.startswith('<Spectrum(flux=<Quantity [')
    assert 'spectral_axis=<SpectralAxis' in result

    spec_with_unc = Spectrum(spectral_axis=wav, flux=flux,
                             uncertainty=StdDevUncertainty(flux))
    result = repr(spec_with_unc)
    assert result.startswith('<Spectrum(flux=<Quantity [')
    assert 'spectral_axis=<SpectralAxis' in result
    assert 'uncertainty=StdDevUncertainty' in result


def test_str():
    wav = np.linspace(100, 1000, 10) * u.nm
    flux = np.ones(10) * u.Jy
    spec = Spectrum(spectral_axis=wav, flux=flux)
    result = str(spec)
    # Sanity check for contents of string representation
    assert result.startswith('Spectrum (length={})'.format(len(spec.flux)))
    lines = result.split('\n')
    flux = spec.flux
    sa = spec.spectral_axis
    assert len(lines) == 3
    assert lines[1].startswith('Flux=')
    assert 'mean={:.5f}'.format(np.nanmean(flux)) in lines[1]
    assert lines[2].startswith('Spectral Axis=')
    assert f'mean={np.nanmean(sa):.5f}' in lines[2]

    # Test string representation with uncertainty
    spec_with_uncertainty = Spectrum(spectral_axis=wav,
                                     flux=flux,
                                     uncertainty=StdDevUncertainty(flux.value))
    result = str(spec_with_uncertainty)
    lines = result.split('\n')
    unc = spec_with_uncertainty.uncertainty
    print(spec_with_uncertainty)
    assert len(lines) == 4
    assert lines[3].startswith('Uncertainty')
    assert f'StdDevUncertainty ([{unc.array[0]:.0f}' in lines[3]

    # Test string representation with multiple flux
    spec_multi_flux = Spectrum(spectral_axis=wav,
                               flux=np.ones((3, 10)) * u.Jy)
    result = str(spec_multi_flux)
    lines = result.split('\n')
    assert len(lines) == 5
    assert lines[1].startswith('Flux=')

    # Test string representation with single-dimensional flux
    spec_single_flux = Spectrum([1] * u.Jy, [0] * u.nm)
    result = str(spec_single_flux)
    assert result == \
"""Spectrum (length=1)
Flux=[1.] Jy,  mean=1.00000 Jy
Spectral Axis=[0.] nm,  mean=0.00000 nm"""


def test_equivalencies():
    """
    Test that after import `u.spectral` equivalencies are not enabled in the global namespace.
    """
    assert u.micron.is_equivalent(u.cm**-1) is False
    assert u.micron.is_equivalent(u.Hz) is False
    assert u.micron.is_equivalent(u.eV) is False
    assert u.Hz.is_equivalent(u.cm**-1) is False
    assert u.Hz.is_equivalent(u.eV) is False


def test_collapse_flux():
    flux = [[2, 4, 6], [0, 8, 12]] * u.Jy
    sa = [100, 200, 300] * u.um
    mask = [[False, True, False], [True, False, False]]
    spec = Spectrum(flux, sa, mask=mask, spectral_axis_index=1)

    assert spec.mean() == 7 * u.Jy
    assert spec.max() == 12 * u.Jy
    assert spec.min() == 2 * u.Jy
    assert spec.sum() == 28 * u.Jy
    assert spec.median() == 7 * u.Jy

    sum_spec = spec.sum(axis = 0)
    assert isinstance(sum_spec, Spectrum)
    assert np.all(sum_spec.flux == [2, 8, 18] * u.Jy)

    mean_spec = spec.mean(axis = 0)
    assert isinstance(mean_spec, Spectrum)
    assert np.all(mean_spec.flux == [2, 8, 9] * u.Jy)

    max_spec = spec.max(axis = 0)
    assert isinstance(max_spec, Spectrum)
    assert np.all(max_spec.flux == [2, 8, 12] * u.Jy)

    min_spec = spec.min(axis = 0)
    assert isinstance(min_spec, Spectrum)
    assert np.all(min_spec.flux == [2, 8, 6] * u.Jy)

    median_spec = spec.mean(axis = 0)
    assert isinstance(median_spec, Spectrum)
    assert np.all(median_spec.flux == [2, 8, 9] * u.Jy)


def test_unsorted_spectral_axis_fails():
    """
    Test that creating a Spectrum fails if spectral axis isn't strictly
    ascending or descending.
    """
    wave = [2, 1, 3] * u.nm
    flux = [5, 5, 5] * u.Jy

    with pytest.raises(ValueError, match='Spectral axis must be strictly increasing or decreasing.'):
        Spectrum(spectral_axis=wave, flux=flux)


def test_spectral_axis_direction():
    """
    Test that the spec1d.spectral_axis_direction attribute correctly reflects
    if the spectral axis in Spectrum is increasing or decreasing.
    """
    flux = [5, 5, 5] * u.Jy

    wave = [1, 2, 3] * u.nm
    spec1d = Spectrum(spectral_axis=wave, flux=flux)
    assert spec1d.spectral_axis_direction == 'increasing'

    wave = [3, 2, 1] * u.nm
    spec1d = Spectrum(spectral_axis=wave, flux=flux)
    assert spec1d.spectral_axis_direction == 'decreasing'
