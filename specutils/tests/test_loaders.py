import logging
import os
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
from astropy.io.registry import IORegistryError
from astropy.modeling import models
from astropy.tests.helper import quantity_allclose

from numpy.testing import assert_allclose

from .conftest import remote_data_path, remote_access
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
def test_spectrumlist_GMOSfits(remote_data_path, caplog):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', (VerifyWarning, UnitsWarning))
        spectrum_list = SpectrumList.read(remote_data_path, format='wcs1d-fits')

    assert len(spectrum_list) == 1

    spec = spectrum_list[0]
    assert len(spec.data) == 3020

    logmsg = caplog.record_tuples[0]
    assert logmsg[1] == logging.WARN
    assert "Assuming the axis 0 labeled 'linear' is spectral" in logmsg[2]



@remote_access([{'id': '1481190', 'filename': 'L5g_0355+11_Cruz09.fits'}])
def test_specific_spec_axis_unit(remote_data_path):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', (VerifyWarning, UnitsWarning))
        optical_spec = Spectrum1D.read(remote_data_path,
                                       spectral_axis_unit="Angstrom",
                                       format='wcs1d-fits')

    assert optical_spec.spectral_axis.unit == "Angstrom"

@remote_access([{'id': '2656720', 'filename': '_v1410ori_20181204_261_Forrest%20Sims.fit'}])
def test_ctypye_not_compliant(remote_data_path, caplog):
    optical_spec = Spectrum1D.read(remote_data_path,
                                   spectral_axis_unit="Angstrom",
                                   format='wcs1d-fits')
    logmsg = caplog.record_tuples[0]
    assert logmsg[1] == logging.WARN
    assert "Assuming the axis 0 labeled 'wavelength' is spectral" in logmsg[2]


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


@pytest.mark.parametrize("name", ['file.fit', 'file.fits', 'file.dat'])
def test_no_reader_matches(name):
    '''If no reader matches a file, check that the correct error is raised.
    This test serves a second purpose: A badly written identifier
    function might raise an error as supposed to returning False when
    it cannot identify a file.  The fact that this test passes means
    that at the very least all identifier functions that have been
    tried for that file ending did not fail with an error.
    '''
    with tempfile.TemporaryDirectory() as tmpdirname:
        filename = os.path.join(tmpdirname, name)
        with open(filename, 'w') as fp:
            fp.write('asdfadasdadvzxcv')

        with pytest.raises(IORegistryError):
            spec = Spectrum1D.read(filename)


@remote_access([{'id':'3359174', 'filename':'linear_fits_solution.fits'}])
def test_iraf_linear(remote_data_path):

    spectrum_1d = Spectrum1D.read(remote_data_path, format='iraf')

    assert isinstance(spectrum_1d, Spectrum1D)
    assert quantity_allclose(spectrum_1d.wavelength[0],
                             u.Quantity(3514.56625402, unit='Angstrom'))
    assert quantity_allclose(spectrum_1d.wavelength[100],
                             u.Quantity(3514.56625402, unit='Angstrom') +
                             u.Quantity(0.653432383823 * 100, unit='Angstrom'))


@remote_access([{'id':'3359180', 'filename':'log-linear_fits_solution.fits'}])
def test_iraf_log_linear(remote_data_path):

    with pytest.raises(NotImplementedError):
        assert Spectrum1D.read(remote_data_path, format='iraf')


@remote_access([{'id':'3359190', 'filename':'non-linear_fits_solution_cheb.fits'}])
def test_iraf_non_linear_chebyshev(remote_data_path):
    chebyshev_model = models.Chebyshev1D(degree=2, domain=[1616, 3259])
    chebyshev_model.c0.value = 5115.64008186
    chebyshev_model.c1.value = 535.515983712
    chebyshev_model.c2.value = -0.779265625182

    wavelength_axis = chebyshev_model(range(1, 4097)) * u.angstrom

    spectrum_1d = Spectrum1D.read(remote_data_path, format='iraf')

    assert isinstance(spectrum_1d, Spectrum1D)
    assert_allclose(wavelength_axis, spectrum_1d.wavelength)


@remote_access([{'id':'3359194', 'filename':'non-linear_fits_solution_legendre.fits'}])
def test_iraf_non_linear_legendre(remote_data_path):

    legendre_model = models.Legendre1D(degree=3, domain=[21, 4048])
    legendre_model.c0.value = 5468.67555891
    legendre_model.c1.value = 835.332144466
    legendre_model.c2.value = -6.02202094803
    legendre_model.c3.value = -1.13142953897

    wavelength_axis = legendre_model(range(1, 4143)) * u.angstrom

    spectrum_1d = Spectrum1D.read(remote_data_path, format='iraf')

    assert isinstance(spectrum_1d, Spectrum1D)
    assert_allclose(wavelength_axis, spectrum_1d.wavelength)


@remote_access([{'id':'3359196', 'filename':'non-linear_fits_solution_linear-spline.fits'}])
def test_iraf_non_linear_linear_spline(remote_data_path):

    with pytest.raises(NotImplementedError):
        assert Spectrum1D.read(remote_data_path, format='iraf')


@remote_access([{'id':'3359200', 'filename':'non-linear_fits_solution_cubic-spline.fits'}])
def test_iraf_non_linear_cubic_spline(remote_data_path):

    with pytest.raises(NotImplementedError):
        assert Spectrum1D.read(remote_data_path, format='iraf')
