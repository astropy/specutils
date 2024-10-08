"""This test module exists because Jdaviz wanted to put pixel area
in the flux unit. Some code copied over from original Jdaviz
implementation.

"""
import numpy as np
import pytest
from astropy import units as u
from astropy.nddata import StdDevUncertainty
from astropy.tests.helper import assert_quantity_allclose
from astropy.wcs import WCS
from numpy.testing import assert_array_equal

from specutils import Spectrum1D, SpectralAxis

PIX2 = u.pix * u.pix


def _eqv_flux_to_sb_pixel():
    """This allows conversion between flux and flux-per-square-pixel
    surface brightness, e.g., MJy and MJy/PIX2.
    """

    # generate an equivalency for each flux type that would need
    # another equivalency for converting to/from
    flux_units = [u.MJy, u.erg / (u.s * u.cm**2 * u.Angstrom),
                  u.ph / (u.Angstrom * u.s * u.cm**2),
                  u.ph / (u.Hz * u.s * u.cm**2)]
    return [(flux_unit, flux_unit / PIX2, lambda x: x, lambda x: x)
            for flux_unit in flux_units]


# The original Jdaviz implementation we are replacing with native
# specutils functionality.
def convert_spectrum1d_from_flux_to_flux_per_pixel(spectrum):
    """Converts a Spectrum1D object's flux units to flux per square pixel.

    This function takes a `specutils.Spectrum1D` object with flux units and converts the
    flux (and optionally, uncertainty) to a surface brightness per square pixel
    (e.g., from Jy to Jy/pix**2). This is done by updating the units of spectrum.flux
    and (if present) spectrum.uncertainty, and creating a new `specutils.Spectrum1D`
    object with the modified flux and uncertainty.

    Parameters
    ----------
    spectrum : Spectrum1D
        A `specutils.Spectrum1D` object containing flux data, which is assumed to be in
        flux units without any angular component in the denominator.

    Returns
    -------
    Spectrum1D
        A new `specutils.Spectrum1D` object with flux and uncertainty (if present)
        converted to units of flux per square pixel.

    """
    # convert flux, which is always populated
    flux = getattr(spectrum, 'flux')
    flux = flux / PIX2

    # and uncerts, if present
    uncerts = getattr(spectrum, 'uncertainty')
    if uncerts is not None:
        # enforce common uncert type.
        uncerts = uncerts.represent_as(StdDevUncertainty)
        uncerts = StdDevUncertainty(uncerts.quantity / PIX2)

    # create a new spectrum 1d with all the info from the input spectrum 1d,
    # and the flux / uncerts converted from flux to SB per square pixel

    # if there is a spectral axis that is a SpectralAxis, you cant also set
    # redshift or radial_velocity
    spectral_axis = getattr(spectrum, 'spectral_axis', None)
    if spectral_axis is not None:
        if isinstance(spectral_axis, SpectralAxis):
            redshift = None
            radial_velocity = None
        else:
            redshift = spectrum.redshift
            radial_velocity = spectrum.radial_velocity

    # initialize new spectrum1d with new flux, uncerts, and all other init parameters
    # from old input spectrum as well as any 'meta'. any more missing information
    # not in init signature that might be present in `spectrum`?
    new_spec1d = Spectrum1D(flux=flux, uncertainty=uncerts,
                            spectral_axis=spectrum.spectral_axis,
                            mask=spectrum.mask,
                            wcs=spectrum.wcs,
                            velocity_convention=spectrum.velocity_convention,
                            rest_value=spectrum.rest_value, redshift=redshift,
                            radial_velocity=radial_velocity,
                            bin_specification=getattr(spectrum, 'bin_specification', None),
                            meta=spectrum.meta)

    return new_spec1d


def assert_dict_equal(d1, d2):
    keys = sorted(d1)
    assert sorted(d2) == keys
    for k in keys:
        v1 = d1[k]
        v2 = d2[k]
        assert v1 == v2


@pytest.mark.parametrize("suppress_conversion", [False, True])
def test_spec_flux_conv_pix2(suppress_conversion):
    meta = {"CTYPE1": "WAVE-LOG", "CTYPE2": "DEC--TAN", "CTYPE3": "RA---TAN",
            "CRVAL1": 4.622e-7, "CRVAL2": 27, "CRVAL3": 205,
            "CDELT1": 8e-11, "CDELT2": 0.0001, "CDELT3": -0.0001,
            "CRPIX1": 0, "CRPIX2": 0, "CRPIX3": 0, "PIXAR_SR": 8e-11}
    w = WCS(meta)
    flux_orig = np.arange(24).reshape((2, 3, 4)) * u.Jy
    uncert_orig = StdDevUncertainty(flux_orig)
    mask_orig = np.zeros(flux_orig.shape, dtype=bool)
    rs_orig = 0.0001 * u.dimensionless_unscaled
    sp_orig = Spectrum1D(
        flux=flux_orig, uncertainty=uncert_orig, mask=mask_orig, wcs=w,
        redshift=rs_orig, meta=meta)

    # What Jdaviz implemented.
    sp_pix2_jdaviz = convert_spectrum1d_from_flux_to_flux_per_pixel(sp_orig)

    # What we should replace it with using only specutils.
    sp_pix2_specutils = sp_orig.with_flux_unit(
        u.Jy / PIX2, equivalencies=_eqv_flux_to_sb_pixel(),
        suppress_conversion=suppress_conversion)

    # Make sure the two implementations are equivalent.
    assert sp_pix2_specutils.flux.unit == u.Jy / PIX2
    assert sp_pix2_specutils.uncertainty.unit == u.Jy / PIX2
    assert_quantity_allclose(sp_pix2_specutils.flux, sp_pix2_jdaviz.flux)
    assert_quantity_allclose(
        sp_pix2_specutils.uncertainty.quantity, sp_pix2_jdaviz.uncertainty.quantity)
    assert_quantity_allclose(sp_pix2_specutils.spectral_axis, sp_pix2_jdaviz.spectral_axis)
    assert_array_equal(sp_pix2_specutils.mask, sp_pix2_jdaviz.mask)
    assert_dict_equal(sp_pix2_specutils.wcs.to_header(), sp_pix2_jdaviz.wcs.to_header())
    assert_dict_equal(sp_pix2_specutils.meta, sp_pix2_jdaviz.meta)
    assert sp_pix2_specutils.velocity_convention == sp_pix2_jdaviz.velocity_convention
    assert sp_pix2_specutils.rest_value is sp_pix2_jdaviz.rest_value  # None
    assert_quantity_allclose(sp_pix2_specutils.redshift, sp_pix2_jdaviz.redshift)
    assert_quantity_allclose(sp_pix2_specutils.radial_velocity, sp_pix2_jdaviz.radial_velocity)
    if "bin_specification" in sp_orig:
        assert sp_pix2_specutils.bin_specification == sp_pix2_jdaviz.bin_specification
