import pickle

import pytest
import numpy as np
from astropy import units as u
from astropy import modeling

from specutils.utils import QuantityModel
from ..utils.wcs_utils import refraction_index, vac_to_air, air_to_vac

wavelengths = [300, 500, 1000] * u.nm
data_index_refraction = {
    'Greisen2006': np.array([3.07393068, 2.9434858 , 2.8925797 ]),
    'Edlen1953': np.array([2.91557413, 2.78963801, 2.74148172]),
    'Edlen1966': np.array([2.91554272, 2.7895973 , 2.74156098]),
    'PeckReeder1972': np.array([2.91554211, 2.78960005, 2.74152561]),
    'Morton2000': np.array([2.91568573, 2.78973402, 2.74169531]),
    'Ciddor1996': np.array([2.91568633, 2.78973811, 2.74166131])
}


def test_quantity_model():
    c = modeling.models.Chebyshev1D(3)
    uc = QuantityModel(c, u.AA, u.km)

    assert uc(10*u.nm).to(u.m) == 0*u.m


def test_pickle_quantity_model(tmp_path):
    """
    Check that a QuantityModel can roundtrip through pickling, as it
    would if fit in a multiprocessing pool.
    """

    c = modeling.models.Chebyshev1D(3)
    uc = QuantityModel(c, u.AA, u.km)

    pkl_file = tmp_path / "qmodel.pkl"

    with open(pkl_file, "wb") as f:
        pickle.dump(uc, f)

    with open(pkl_file, "rb") as f:
        new_model = pickle.load(f)

    assert new_model.input_units == uc.input_units
    assert new_model.return_units == uc.return_units
    assert isinstance(new_model.unitless_model, uc.unitless_model.__class__)
    assert np.all(new_model.unitless_model.parameters == uc.unitless_model.parameters)


@pytest.mark.parametrize("method", data_index_refraction.keys())
def test_refraction_index(method):
    tmp = (refraction_index(wavelengths, method) - 1) * 1e4
    assert np.isclose(tmp, data_index_refraction[method], atol=1e-7).all()


@pytest.mark.parametrize("method", data_index_refraction.keys())
def test_air_to_vac(method):
    tmp = refraction_index(wavelengths, method)
    assert np.isclose(wavelengths.value * tmp,
                      air_to_vac(wavelengths, method=method, scheme='inversion').value,
                      rtol=1e-6).all()
    assert np.isclose(wavelengths.value,
                      air_to_vac(vac_to_air(wavelengths, method=method),
                                 method=method, scheme='iteration').value,
                      atol=1e-12).all()


def test_refraction_index_conditions():
    standard = refraction_index(wavelengths)

    # Standard conditions leave the result unchanged, in either temperature unit
    assert np.allclose(refraction_index(wavelengths, temperature=15 * u.deg_C,
                                        pressure=101325 * u.Pa), standard, rtol=0, atol=1e-15)
    assert np.allclose(refraction_index(wavelengths, temperature=288.15 * u.K),
                       standard, rtol=0, atol=1e-15)

    # n - 1 scales roughly with density: lower pressure and higher temperature
    # both reduce it
    low_pressure = refraction_index(wavelengths, pressure=700 * u.hPa, temperature=5 * u.deg_C)
    density_ratio = (700e2 / 101325) * (288.15 / 278.15)
    assert np.allclose((low_pressure - 1) / (standard - 1), density_ratio, rtol=2e-3)
    assert np.all(refraction_index(wavelengths, temperature=30 * u.deg_C) < standard)

    # Water vapour lowers the refractive index slightly, by ~3e-7 at 50% humidity
    humid = refraction_index(wavelengths, humidity=0.5)
    assert np.all(humid < standard)
    assert np.allclose(standard - humid, 3e-7, rtol=0.1)
    assert np.allclose(refraction_index(wavelengths, humidity=50 * u.percent), humid)

    with pytest.raises(ValueError, match="humidity must be"):
        refraction_index(wavelengths, humidity=1.5)
    with pytest.raises(u.UnitConversionError):
        refraction_index(wavelengths, pressure=5 * u.K)


def test_air_to_vac_conditions():
    conditions = dict(temperature=5 * u.deg_C, pressure=700 * u.hPa, humidity=0.2)
    air = vac_to_air(wavelengths, **conditions)
    assert np.all(air > vac_to_air(wavelengths))  # thinner air refracts less
    for scheme in ('inversion', 'iteration'):
        assert np.allclose(air_to_vac(air, scheme=scheme, **conditions).value,
                           wavelengths.value, atol=1e-6)
    with pytest.raises(ValueError, match="not supported with scheme='Piskunov'"):
        air_to_vac(air, scheme='Piskunov', **conditions)
