import shutil
import tempfile
import urllib
import warnings

import pytest

import astropy.units as u
import numpy as np
from astropy.io.fits.verify import VerifyWarning
from astropy.table import Table
from astropy.units import UnitsWarning
from astropy.wcs import FITSFixedWarning

from .conftest import remote_access
from .. import Spectrum1D, SpectrumList
from ..io import get_loaders_by_extension


def test_get_loaders_by_extension():
    loader_labels = get_loaders_by_extension('fits')

    assert len(loader_labels) > 0
    assert isinstance(loader_labels[0], str)


@remote_access([{'id': '1481190', 'filename': 'L5g_0355+11_Cruz09.fits'}])
def test_spectrum1d_GMOSfits(remote_data_path):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', (VerifyWarning, UnitsWarning))
        optical_spec_2 = Spectrum1D.read(remote_data_path, format='wcs1d-fits')

    assert len(optical_spec_2.data) == 3020


@remote_access([{'id': '1481190', 'filename': 'L5g_0355+11_Cruz09.fits'}])
def test_spectrumlist_GMOSfits(remote_data_path):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', (VerifyWarning, UnitsWarning))
        spectrum_list = SpectrumList.read(remote_data_path, format='wcs1d-fits')

    assert len(spectrum_list) == 1

    spec = spectrum_list[0]
    assert len(spec.data) == 3020


@remote_access([{'id': '1481190', 'filename': 'L5g_0355+11_Cruz09.fits'}])
def test_specific_spec_axis_unit(remote_data_path):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', (VerifyWarning, UnitsWarning))
        optical_spec = Spectrum1D.read(remote_data_path,
                                       spectral_axis_unit="Angstrom",
                                       format='wcs1d-fits')

    assert optical_spec.spectral_axis.unit == "Angstrom"


def test_generic_ecsv_reader(tmpdir):
   # Create a small data set
   wave = np.arange(1,1.1,0.01)*u.AA
   flux = np.ones(len(wave))*1.e-14*u.Jy
   uncertainty = 0.01*flux
   table = Table([wave,flux,uncertainty],names=["wave","flux","uncertainty"])
   tmpfile = str(tmpdir.join('_tst.ecsv'))
   table.write(tmpfile,format='ascii.ecsv')

   # Read it in and check against the original
   spectrum = Spectrum1D.read(tmpfile,format='ECSV')
   assert spectrum.spectral_axis.unit == table['wave'].unit
   assert spectrum.flux.unit == table['flux'].unit
   assert spectrum.uncertainty.unit == table['uncertainty'].unit
   assert spectrum.spectral_axis.unit == table['wave'].unit
   assert np.alltrue(spectrum.spectral_axis == table['wave'])
   assert np.alltrue(spectrum.flux == table['flux'])
   assert np.alltrue(spectrum.uncertainty.array == table['uncertainty'])


@remote_access([{'id': '1481119', 'filename': 'COS_FUV.fits'},
                {'id': '1481181', 'filename': 'COS_NUV.fits'}])
def test_hst_cos(remote_data_path):
    spec = Spectrum1D.read(remote_data_path, format='HST/COS')

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0


@remote_access([{'id': '1481192', 'filename':'STIS_FUV.fits'},
                {'id': '1481185', 'filename': 'STIS_NUV.fits'},
                {'id': '1481183', 'filename': 'STIS_CCD.fits'}])
def test_hst_stis(remote_data_path):
    spec = Spectrum1D.read(remote_data_path, format='HST/STIS')

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0


@pytest.mark.remote_data
def test_sdss_spec():
    with urllib.request.urlopen('https://dr14.sdss.org/optical/spectrum/view/data/format%3Dfits/spec%3Dlite?mjd=55359&fiberid=596&plateid=4055') as response:
        with tempfile.NamedTemporaryFile() as tmp_file:
            shutil.copyfileobj(response, tmp_file)

            spec = Spectrum1D.read(tmp_file.name, format="SDSS-III/IV spec")

            assert isinstance(spec, Spectrum1D)
            assert spec.flux.size > 0


@pytest.mark.remote_data
def test_sdss_spspec():
    with urllib.request.urlopen('http://das.sdss.org/spectro/1d_26/0273/1d/spSpec-51957-0273-016.fit') as response:
        with tempfile.NamedTemporaryFile() as tmp_file:
            shutil.copyfileobj(response, tmp_file)

            with warnings.catch_warnings():
                warnings.simplefilter('ignore', FITSFixedWarning)
                spec = Spectrum1D.read(tmp_file.name, format="SDSS-I/II spSpec")

            assert isinstance(spec, Spectrum1D)
            assert spec.flux.size > 0
