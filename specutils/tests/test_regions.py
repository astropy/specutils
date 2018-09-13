import numpy as np

import astropy.units as u
from astropy.nddata import StdDevUncertainty

from ..spectra import Spectrum1D, SpectralRegion


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
    sr = SpectralRegion([(0.3*u.um, 1.0*u.um), (0.45*u.um, 0.6*u.um), (0.04*u.um, 0.05*u.um),
                         (0.8*u.um, 0.9*u.um)])

    assert sr.lower == 0.04*u.um
    assert sr.upper == 1.0*u.um

    # Get lower bound of a single sub-region:
    assert sr[0].lower == 0.04*u.um
    assert sr[0].upper == 0.05*u.um


def test_adding_spectral_regions():

    # Combine two Spectral regions into one:
    sr = SpectralRegion(0.45*u.um, 0.6*u.um) + SpectralRegion(0.8*u.um, 0.9*u.um)

    assert set(sr.subregions) == set([(0.45*u.um, 0.6*u.um), (0.8*u.um, 0.9*u.um)])

    # In-place adding spectral regions:
    sr1 = SpectralRegion(0.45*u.um, 0.6*u.um)
    sr2 = SpectralRegion(0.8*u.um, 0.9*u.um)
    sr1 += sr2

    assert set(sr1.subregions) == set([(0.45*u.um, 0.6*u.um), (0.8*u.um, 0.9*u.um)])


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

    sr = SpectralRegion(0.15*u.um, 0.2*u.um) + SpectralRegion(0.3*u.um, 0.4*u.um) +\
         SpectralRegion(0.45*u.um, 0.6*u.um) + SpectralRegion(0.8*u.um, 0.9*u.um) +\
         SpectralRegion(1.0*u.um, 1.2*u.um) + SpectralRegion(1.3*u.um, 1.5*u.um)

    subsr = sr[3:5]

    assert subsr[0].subregions == [(0.8*u.um, 0.9*u.um)]
    assert subsr[1].subregions == [(1.0*u.um, 1.2*u.um)]


def test_invert():
    sr = SpectralRegion(0.15*u.um, 0.2*u.um) + SpectralRegion(0.3*u.um, 0.4*u.um) +\
         SpectralRegion(0.45*u.um, 0.6*u.um) + SpectralRegion(0.8*u.um, 0.9*u.um) +\
         SpectralRegion(1.0*u.um, 1.2*u.um) + SpectralRegion(1.3*u.um, 1.5*u.um)

    sr_inverted_expected = [(0.05*u.um, 0.15*u.um), (0.2*u.um, 0.3*u.um), (0.4*u.um, 0.45*u.um),
                            (0.6*u.um, 0.8*u.um), (0.9*u.um, 1.0*u.um), (1.2*u.um, 1.3*u.um),
                            (1.5*u.um, 3.0*u.um)]

    # Invert from range.
    sr_inverted = sr.invert(0.05*u.um, 3*u.um)

    for ii, expected in enumerate(sr_inverted_expected):
        assert sr_inverted.subregions[ii] == sr_inverted_expected[ii]

    # Invert from spectrum.
    spectrum = Spectrum1D(spectral_axis=np.linspace(0.05, 3, 20)*u.um, flux=np.random.random(20)*u.Jy)
    sr_inverted = sr.invert_from_spectrum(spectrum)
    for ii, expected in enumerate(sr_inverted_expected):
        assert sr_inverted.subregions[ii] == sr_inverted_expected[ii]
