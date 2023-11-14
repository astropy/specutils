# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests SpecUtils io routines
"""
from collections import Counter

import astropy.units as u
import numpy as np
import pytest
from astropy.io import registry
from astropy.table import Table
from astropy.utils.exceptions import AstropyUserWarning

from specutils import Spectrum1D, SpectrumList
from specutils.io import data_loader
from specutils.io.parsing_utils import generic_spectrum_from_table  # or something like that
from specutils.io.registers import _astropy_has_priorities


def test_generic_spectrum_from_table(recwarn):
    """
    Read a simple table with wavelength, flux and uncertainty
    """
    # Create a small data set, first without uncertainties
    wave = np.arange(1, 1.1, 0.01)*u.AA
    flux = np.ones(len(wave))*1.e-14*u.Jy
    table = Table([wave, flux], names=["wave", "flux"])

    # Test that the units and values of the Spectrum1D object match those in the table
    spectrum = generic_spectrum_from_table(table)
    assert spectrum.spectral_axis.unit == table['wave'].unit
    assert spectrum.flux.unit == table['flux'].unit
    assert spectrum.spectral_axis.unit == table['wave'].unit
    assert np.all(spectrum.spectral_axis == table['wave'])
    assert np.all(spectrum.flux == table['flux'])

    # Add uncertainties and retest
    err = 0.01*flux
    table = Table([wave, flux, err], names=["wave", "flux", "err"])
    spectrum = generic_spectrum_from_table(table)
    assert spectrum.spectral_axis.unit == table['wave'].unit
    assert spectrum.flux.unit == table['flux'].unit
    assert spectrum.uncertainty.unit == table['err'].unit
    assert spectrum.spectral_axis.unit == table['wave'].unit
    assert np.all(spectrum.spectral_axis == table['wave'])
    assert np.all(spectrum.flux == table['flux'])
    assert np.all(spectrum.uncertainty.array == table['err'])

    # Test for warning if standard deviation is zero or negative
    err[0] = 0.
    table = Table([wave, flux, err], names=["wave", "flux", "err"])
    spectrum = generic_spectrum_from_table(table)
    assert len(recwarn) == 1
    w = recwarn.pop(AstropyUserWarning)
    assert "Standard Deviation has values of 0 or less" in str(w.message)

    # Test that exceptions are raised if there are no units
    flux = np.ones(len(wave))*1.e-14
    table = Table([wave, flux], names=["wave", "flux"])
    with pytest.raises(IOError) as exc:
        spectrum = generic_spectrum_from_table(table)
        assert 'Could not identify column containing the flux' in exc
    wave = np.arange(1, 1.1, 0.01)
    table = Table([wave, flux, err], names=["wave", "flux", "err"])
    with pytest.raises(IOError) as exc:
        spectrum = generic_spectrum_from_table(table)
        assert 'Could not identify column containing the wavelength, frequency or energy' in exc


def test_speclist_autoidentify():

    formats = registry.get_formats(SpectrumList)
    assert (formats['Auto-identify'] == 'Yes').all()


@pytest.mark.filterwarnings(r'ignore:.*data loader provided for Spectrum1D without explicit identifier')
def test_default_identifier(tmp_path):

    fname = str(tmp_path / 'empty.txt')
    with open(fname, 'w') as ff:
        ff.write('\n')

    format_name = 'default_identifier_test'

    @data_loader(format_name)
    def reader(*args, **kwargs):
        """Doesn't actually get used."""
        return

    for datatype in [Spectrum1D, SpectrumList]:
        fmts = registry.identify_format('read', datatype, fname, None, [], {})
        assert format_name in fmts

        # Clean up after ourselves
        registry.unregister_reader(format_name, datatype)
        registry.unregister_identifier(format_name, datatype)


def test_default_identifier_extension(tmp_path):

    good_fname = str(tmp_path / 'empty.fits')
    bad_fname = str(tmp_path / 'empty.txt')

    # Create test data files.
    for name in [good_fname, bad_fname]:
        with open(name, 'w') as ff:
            ff.write('\n')

    format_name = 'default_identifier_extension_test'

    @data_loader(format_name, extensions=['fits'])
    def reader(*args, **kwargs):
        """Doesn't actually get used."""
        return

    for datatype in [Spectrum1D, SpectrumList]:
        fmts = registry.identify_format('read', datatype, good_fname, None, [], {})
        assert format_name in fmts

        fmts = registry.identify_format('read', datatype, bad_fname, None, [], {})
        assert format_name not in fmts

        # Clean up after ourselves
        registry.unregister_reader(format_name, datatype)
        registry.unregister_identifier(format_name, datatype)


def test_custom_identifier(tmp_path):

    good_fname = str(tmp_path / 'good.txt')
    bad_fname = str(tmp_path / 'bad.txt')

    # Create test data files.
    for name in [good_fname, bad_fname]:
        with open(name, 'w') as ff:
            ff.write('\n')

    format_name = 'custom_identifier_test'

    def identifier(origin, *args, **kwargs):
        fname = args[0]
        return 'good' in fname

    @data_loader(format_name, identifier=identifier)
    def reader(*args, **kwargs):
        """Doesn't actually get used."""
        return

    for datatype in [Spectrum1D, SpectrumList]:
        fmts = registry.identify_format('read', datatype, good_fname, None, [], {})
        assert format_name in fmts

        fmts = registry.identify_format('read', datatype, bad_fname, None, [], {})
        assert format_name not in fmts

        # Clean up after ourselves
        registry.unregister_reader(format_name, datatype)
        registry.unregister_identifier(format_name, datatype)


@pytest.mark.xfail(
    not _astropy_has_priorities(),
    reason="Test requires priorities to be implemented in astropy",
    raises=registry.IORegistryError,
)
def test_loader_uses_priority(tmp_path):
    counter = Counter()
    fname = str(tmp_path / 'good.txt')

    with open(fname, 'w') as ff:
        ff.write('\n')

    def identifier(origin, *args, **kwargs):
        fname = args[0]
        return 'good' in fname

    @data_loader("test_counting_loader1", identifier=identifier, priority=1)
    def counting_loader1(*args, **kwargs):
        counter["test1"] += 1
        wave = np.arange(1, 1.1, 0.01)*u.AA
        return Spectrum1D(
            spectral_axis=wave,
            flux=np.ones(len(wave))*1.e-14*u.Jy,
        )

    @data_loader("test_counting_loader2", identifier=identifier, priority=2)
    def counting_loader2(*args, **kwargs):
        counter["test2"] += 1
        wave = np.arange(1, 1.1, 0.01)*u.AA
        return Spectrum1D(
            spectral_axis=wave,
            flux=np.ones(len(wave))*1.e-14*u.Jy,
        )

    Spectrum1D.read(fname)
    assert counter["test2"] == 1
    assert counter["test1"] == 0

    for datatype in [Spectrum1D, SpectrumList]:
        registry.unregister_reader("test_counting_loader1", datatype)
        registry.unregister_identifier("test_counting_loader1", datatype)
        registry.unregister_reader("test_counting_loader2", datatype)
        registry.unregister_identifier("test_counting_loader2", datatype)
