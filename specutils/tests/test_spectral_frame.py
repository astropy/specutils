import astropy.units as u
import pytest

from ..spectra.spectral_frame import SpectralMedium, SPECTRAL_FRAMES, normalize_frame


def test_spectral_medium_vacuum():
    medium = SpectralMedium('vacuum')
    assert medium.is_vacuum and not medium.is_air
    assert medium.refraction_method is None
    assert medium.refraction_kwargs == {}
    assert str(medium) == 'vacuum'
    assert repr(medium) == "SpectralMedium(kind='vacuum')"
    assert medium == SpectralMedium('VACUUM')

    with pytest.raises(ValueError, match="Air conditions cannot be set"):
        SpectralMedium('vacuum', co2=400)


def test_spectral_medium_air():
    medium = SpectralMedium('air')
    assert medium.is_air
    assert medium.refraction_method == 'Morton2000'
    assert medium.refraction_kwargs['method'] == 'Morton2000'

    medium = SpectralMedium('air', refraction_method='Ciddor1996', co2=400,
                            temperature=10 * u.deg_C, pressure=700 * u.hPa,
                            humidity=20 * u.percent)
    assert medium.co2 == 400
    assert medium.temperature == 10 * u.deg_C
    assert medium.pressure == 700 * u.hPa
    assert medium.humidity == pytest.approx(0.2)
    assert medium.to_dict()['refraction_method'] == 'Ciddor1996'
    assert medium != SpectralMedium('air')


def test_spectral_medium_validation():
    with pytest.raises(ValueError, match="kind must be one of"):
        SpectralMedium('water')
    with pytest.raises(u.UnitsError):
        SpectralMedium('air', temperature=10 * u.m)
    with pytest.raises(u.UnitsError):
        SpectralMedium('air', pressure=10 * u.K)
    with pytest.raises(ValueError, match="humidity must be between"):
        SpectralMedium('air', humidity=1.5)


def test_spectral_medium_from_input():
    assert SpectralMedium.from_input(None) is None
    assert SpectralMedium.from_input('air') == SpectralMedium('air')
    medium = SpectralMedium('air', co2=400)
    assert SpectralMedium.from_input(medium) is medium
    assert SpectralMedium.from_input({'kind': 'air', 'co2': 400}) == medium
    with pytest.raises(TypeError):
        SpectralMedium.from_input(3)


@pytest.mark.parametrize('frame, expected', [
    (None, None),
    ('BARYCENT', 'BARYCENT'),
    ('barycent', 'BARYCENT'),
    ('barycentric', 'BARYCENT'),
    ('topocentric', 'TOPOCENT'),
    ('rest', 'SOURCE'),
    ('source', 'SOURCE'),
    ('lsrk', 'LSRK'),
])
def test_normalize_frame(frame, expected):
    assert normalize_frame(frame) == expected


def test_normalize_frame_invalid():
    with pytest.raises(ValueError, match="Unknown spectral frame"):
        normalize_frame('ecliptic')
    with pytest.raises(TypeError):
        normalize_frame(3)
    for frame in SPECTRAL_FRAMES:
        assert normalize_frame(frame) == frame
