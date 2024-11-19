import astropy.units as u
import numpy as np
import pytest
from astropy.modeling.models import Gaussian1D
from astropy.tests.helper import assert_quantity_allclose

from specutils.fitting import find_lines_derivative
from specutils.spectra import Spectrum1D, SpectralRegion


def test_lower_upper():
    # Spectral region with just one range (lower and upper bound)
    sr = SpectralRegion(0.45*u.um, 0.6*u.um)

    assert sr.lower == 0.45*u.um
    assert sr.upper == 0.6*u.um

    # Spectral region with just two ranges
    sr = SpectralRegion([(0.45*u.um, 0.6*u.um), (0.8*u.um, 0.9*u.um)])

    assert sr.lower == 0.45*u.um
    assert sr.upper == 0.9*u.um

    # Spectral region with multiple ranges and not ordered
    sr = SpectralRegion([(0.3*u.um, 1.0*u.um), (0.45*u.um, 0.6*u.um),
                         (0.04*u.um, 0.05*u.um), (0.8*u.um, 0.9*u.um)])

    assert sr.lower == 0.04*u.um
    assert sr.upper == 1.0*u.um

    # Get lower bound of a single sub-region:
    assert sr[0].lower == 0.04*u.um
    assert sr[0].upper == 0.05*u.um


@pytest.mark.parametrize(
    ('center', 'width', 'lower', 'upper'),
    [(6563 * u.AA, 10 * u.AA, 6558.0 * u.AA, 6568.0 * u.AA),
     (1 * u.GHz, 0.1 * u.GHz, 1.05 * u.GHz, 0.95 * u.GHz),
     (0.5 * u.pix, 1 * u.pix, 0 * u.pix, 1 * u.pix)])
def test_from_center(center, width, lower, upper):
    # Spectral region from center with width
    sr = SpectralRegion.from_center(center=center, width=width)
    assert_quantity_allclose(sr.lower, lower)
    assert_quantity_allclose(sr.upper, upper)


@pytest.mark.parametrize(
    ('center', 'width'),
    [(6563 * u.AA, -10 * u.AA),
     (6563 * u.AA, 0 * u.AA),
     (1 * u.GHz, -0.1 * u.GHz)])
def test_from_center_error(center, width):
    with pytest.raises(ValueError):
        SpectralRegion.from_center(center=center, width=width)


def test_adding_spectral_regions():

    # Combine two Spectral regions into one:
    sr = (SpectralRegion(0.45*u.um, 0.6*u.um) +
          SpectralRegion(0.8*u.um, 0.9*u.um))

    assert set(sr.subregions) == set([(0.45*u.um, 0.6*u.um),
                                      (0.8*u.um, 0.9*u.um)])

    # In-place adding spectral regions:
    sr1 = SpectralRegion(0.45*u.um, 0.6*u.um)
    sr2 = SpectralRegion(0.8*u.um, 0.9*u.um)
    sr1 += sr2

    assert set(sr1.subregions) == set([(0.45*u.um, 0.6*u.um),
                                       (0.8*u.um, 0.9*u.um)])


def test_getitem():
    sr = SpectralRegion([(0.8*u.um, 0.9*u.um), (0.3*u.um, 1.0*u.um),
                         (0.45*u.um, 0.6*u.um), (0.04*u.um, 0.05*u.um)])

    assert sr[0].subregions == [(0.04*u.um, 0.05*u.um)]
    assert sr[1].subregions == [(0.3*u.um, 1.0*u.um)]
    assert sr[2].subregions == [(0.45*u.um, 0.6*u.um)]
    assert sr[3].subregions == [(0.8*u.um, 0.9*u.um)]
    assert sr[-1].subregions == [(0.8*u.um, 0.9*u.um)]


def test_bounds():

    # Single subregion
    sr = SpectralRegion(0.45*u.um, 0.6*u.um)
    assert sr.bounds == (0.45*u.um, 0.6*u.um)

    # Multiple subregions
    sr = SpectralRegion([(0.8*u.um, 0.9*u.um), (0.3*u.um, 1.0*u.um),
                         (0.45*u.um, 0.6*u.um), (0.04*u.um, 0.05*u.um)])
    assert sr.bounds == (0.04*u.um, 1.0*u.um)


def test_delitem():

    # Single subregion
    sr = SpectralRegion(0.45*u.um, 0.6*u.um)

    del sr[0]
    assert sr.subregions == []

    # Multiple sub-regions
    sr = SpectralRegion([(0.8*u.um, 0.9*u.um), (0.3*u.um, 1.0*u.um),
                         (0.45*u.um, 0.6*u.um), (0.04*u.um, 0.05*u.um)])
    del sr[1]

    assert sr[0].subregions == [(0.04*u.um, 0.05*u.um)]
    assert sr[1].subregions == [(0.45*u.um, 0.6*u.um)]
    assert sr[2].subregions == [(0.8*u.um, 0.9*u.um)]


def test_iterate():

    # Create the Spectral region
    subregions = [(0.8*u.um, 0.9*u.um), (0.3*u.um, 1.0*u.um),
                  (0.45*u.um, 0.6*u.um), (0.04*u.um, 0.05*u.um)]
    sr = SpectralRegion(subregions)

    # For testing, sort our subregion list.
    subregions.sort(key=lambda k: k[0])

    for ii, s in enumerate(sr):
        assert s.subregions[0] == subregions[ii]


def test_slicing():

    sr = (SpectralRegion(0.15*u.um, 0.2*u.um) +
          SpectralRegion(0.3*u.um, 0.4*u.um) +
          SpectralRegion(0.45*u.um, 0.6*u.um) +
          SpectralRegion(0.8*u.um, 0.9*u.um) +
          SpectralRegion(1.0*u.um, 1.2*u.um) +
          SpectralRegion(1.3*u.um, 1.5*u.um))

    subsr = sr[3:5]

    assert subsr[0].subregions == [(0.8*u.um, 0.9*u.um)]
    assert subsr[1].subregions == [(1.0*u.um, 1.2*u.um)]


def test_invert():
    sr = (SpectralRegion(0.15*u.um, 0.2*u.um) +
          SpectralRegion(0.3*u.um, 0.4*u.um) +
          SpectralRegion(0.45*u.um, 0.6*u.um) +
          SpectralRegion(0.8*u.um, 0.9*u.um) +
          SpectralRegion(1.0*u.um, 1.2*u.um) +
          SpectralRegion(1.3*u.um, 1.5*u.um))

    sr_inverted_expected = [(0.05*u.um, 0.15*u.um), (0.2*u.um, 0.3*u.um),
                            (0.4*u.um, 0.45*u.um), (0.6*u.um, 0.8*u.um),
                            (0.9*u.um, 1.0*u.um), (1.2*u.um, 1.3*u.um),
                            (1.5*u.um, 3.0*u.um)]

    # Invert from range.
    sr_inverted = sr.invert(0.05*u.um, 3*u.um)

    for ii, expected in enumerate(sr_inverted_expected):
        assert sr_inverted.subregions[ii] == sr_inverted_expected[ii]

    # Invert from spectrum.
    spectrum = Spectrum1D(spectral_axis=np.linspace(0.05, 3, 20)*u.um,
                          flux=np.random.random(20)*u.Jy)
    sr_inverted = sr.invert_from_spectrum(spectrum)
    for ii, expected in enumerate(sr_inverted_expected):
        assert sr_inverted.subregions[ii] == sr_inverted_expected[ii]


def test_from_list_list():
    g1 = Gaussian1D(1, 4.6, 0.2)
    g2 = Gaussian1D(2.5, 5.5, 0.1)
    g3 = Gaussian1D(-1.7, 8.2, 0.1)

    x = np.linspace(0, 10, 200)
    y = g1(x) + g2(x) + g3(x)

    spectrum = Spectrum1D(flux=y * u.Jy, spectral_axis=x * u.um)

    lines = find_lines_derivative(spectrum, flux_threshold=0.01)

    spec_reg = SpectralRegion.from_line_list(lines)
    expected = [(4.072864321608041 * u.um, 5.072864321608041 * u.um),
                (4.977386934673367 * u.um, 5.977386934673367 * u.um),
                (7.690954773869347 * u.um, 8.690954773869347 * u.um)]

    for i, reg in enumerate(expected):
        assert_quantity_allclose(reg, (spec_reg[i].lower, spec_reg[i].upper))


def test_read_write(tmp_path):
    path = tmp_path / "test_sr.ecsv"
    sr = SpectralRegion([(0.45*u.um, 0.6*u.um), (0.8*u.um, 0.9*u.um)])
    sr.write(str(path))
    sr2 = SpectralRegion.read(path)
    assert list(sr2.as_table().columns) == ["lower_bound", "upper_bound"]

    sr3 = SpectralRegion.from_qtable(sr2.as_table())
    assert sr3.subregions[0] == sr.subregions[0]
    assert sr3.subregions[1] == sr.subregions[1]
