import numpy as np
import pytest

import astropy.units as u
from astropy.io import fits
from astropy.nddata import InverseVariance, StdDevUncertainty, VarianceUncertainty
from astropy.tests.helper import quantity_allclose

from specutils import Spectrum


@pytest.mark.parametrize(
    "uncertainty",
    [StdDevUncertainty, VarianceUncertainty, InverseVariance],
)
def test_generic_ecsv_writer(tmp_path, uncertainty):
    spectral_axis = np.linspace(4000, 5000, 11) * u.AA
    flux = np.linspace(1, 2, spectral_axis.size) * u.Jy
    stddev = np.linspace(100, 200, spectral_axis.size) * u.mJy

    if uncertainty is StdDevUncertainty:
        unc = uncertainty(stddev)
    elif uncertainty is VarianceUncertainty:
        unc = uncertainty(stddev**2)
    else:
        unc = uncertainty(stddev**-2)

    mask = np.zeros(flux.shape, dtype=bool)
    mask[3] = True
    header = fits.Header({"OBJECT": "test target"})
    spectrum = Spectrum(
        spectral_axis=spectral_axis,
        flux=flux,
        uncertainty=unc,
        mask=mask,
        meta={"header": header},
    )

    path = tmp_path / "spectrum.ecsv"
    spectrum.write(path)

    result = Spectrum.read(path)
    assert quantity_allclose(result.spectral_axis, spectrum.spectral_axis)
    assert quantity_allclose(result.flux, spectrum.flux)
    expected_uncertainty = (
        spectrum.uncertainty
        .represent_as(StdDevUncertainty)
        .quantity
        .to(flux.unit)
    )
    assert quantity_allclose(result.uncertainty.quantity, expected_uncertainty)
    assert np.array_equal(result.mask, spectrum.mask)
    assert result.meta["header"]["OBJECT"] == "test target"

    with pytest.raises(OSError, match="already exists"):
        spectrum.write(path, format="ECSV")
    spectrum.write(path, format="ECSV", overwrite=True)


def test_generic_ecsv_writer_multidimensional(tmp_path):
    spectral_axis = np.linspace(4000, 5000, 5) * u.AA
    shape = (2, 3, spectral_axis.size)
    flux = np.arange(np.prod(shape)).reshape(shape) * u.Jy
    uncertainty = StdDevUncertainty(np.full(shape, 0.1) * u.Jy)
    mask = np.zeros(shape, dtype=np.uint16)
    mask[0, 1, 2] = 1

    spectrum = Spectrum(
        spectral_axis=spectral_axis,
        flux=flux,
        uncertainty=uncertainty,
        mask=mask,
    )
    path = tmp_path / "multidimensional.ecsv"
    spectrum.write(path)

    result = Spectrum.read(path)
    assert result.flux.shape == spectrum.flux.shape
    assert result.uncertainty.array.shape == spectrum.uncertainty.array.shape
    assert result.mask.shape == spectrum.mask.shape
    assert quantity_allclose(result.spectral_axis, spectrum.spectral_axis)
    assert quantity_allclose(result.flux, spectrum.flux)
    assert quantity_allclose(result.uncertainty.quantity,
                             spectrum.uncertainty.quantity)
    assert np.array_equal(result.mask, spectrum.mask)
