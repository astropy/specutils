import logging
import os
import sys
import shutil
import tempfile
import urllib
import warnings

import pytest

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.table import Table
from astropy.units import UnitsWarning
from astropy.wcs import FITSFixedWarning
from astropy.io.registry import IORegistryError
from astropy.modeling import models
from astropy.tests.helper import quantity_allclose
from astropy.nddata import NDUncertainty, StdDevUncertainty

from numpy.testing import assert_allclose

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
def test_spectrumlist_GMOSfits(remote_data_path, caplog):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', (VerifyWarning, UnitsWarning))
        spectrum_list = SpectrumList.read(remote_data_path,
                                          format='wcs1d-fits')

    assert len(spectrum_list) == 1

    spec = spectrum_list[0]
    assert len(spec.data) == 3020

    assert len(caplog.record_tuples) == 0


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

    assert len(caplog.record_tuples) == 0


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
    spec = Spectrum1D.read(remote_data_path)

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0


@remote_access([{'id': '1481192', 'filename':'STIS_FUV.fits'},
                {'id': '1481185', 'filename': 'STIS_NUV.fits'},
                {'id': '1481183', 'filename': 'STIS_CCD.fits'}])
def test_hst_stis(remote_data_path):
    spec = Spectrum1D.read(remote_data_path)

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0


@pytest.mark.remote_data
def test_sdss_spec():
    sp_pattern = 'spec-4055-55359-0596.fits.'
    with urllib.request.urlopen('https://dr14.sdss.org/optical/spectrum/view/data/format%3Dfits/spec%3Dlite?mjd=55359&fiberid=596&plateid=4055') as response:
        with tempfile.NamedTemporaryFile(prefix=sp_pattern) as tmp_file:
            shutil.copyfileobj(response, tmp_file)

            spec = Spectrum1D.read(tmp_file.name)

            assert isinstance(spec, Spectrum1D)
            assert spec.flux.size > 0


@pytest.mark.remote_data
def test_sdss_spspec():
    sp_pattern = 'spSpec-51957-0273-016.fit.'
    with urllib.request.urlopen('http://das.sdss.org/spectro/1d_26/0273/1d/spSpec-51957-0273-016.fit') as response:
        with tempfile.NamedTemporaryFile(prefix=sp_pattern) as tmp_file:
            shutil.copyfileobj(response, tmp_file)

            with warnings.catch_warnings():
                warnings.simplefilter('ignore', FITSFixedWarning)
                spec = Spectrum1D.read(tmp_file.name, format="SDSS-I/II spSpec")

            assert isinstance(spec, Spectrum1D)
            assert spec.flux.size > 0


@pytest.mark.remote_data
def test_sdss_spec_stream():
    '''Test direct read and recognition of SDSS-III/IV spec from remote URL,
    i.e. do not rely on filename pattern.
    '''
    sdss_url = 'https://dr14.sdss.org/optical/spectrum/view/data/format%3Dfits/spec%3Dlite?mjd=55359&fiberid=596&plateid=4055'
    spec = Spectrum1D.read(sdss_url)

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0
    assert spec.uncertainty.array.min() >= 0.0


@pytest.mark.remote_data
def test_sdss_spspec_stream():
    '''Test direct read and recognition of SDSS-I/II spSpec from remote URL,
    i.e. do not rely on filename pattern.
    '''
    sdss_url = 'http://das.sdss.org/spectro/1d_26/0273/1d/spSpec-51957-0273-016.fit'
    spec = Spectrum1D.read(sdss_url)

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0
    assert spec.uncertainty.array.min() >= 0.0


@pytest.mark.skipif('sys.platform.startswith("win")',
                    reason='Uncertain availability of compression utilities')
@pytest.mark.remote_data
@pytest.mark.parametrize('compress', ['gzip', 'bzip2'])
def test_sdss_compressed(compress):
    '''Test automatic recognition of supported compression formats.
    '''
    ext = {'gzip': '.gz', 'bzip2': '.bz2'}
    # Deliberately not using standard filename pattern to test header info.
    sp_pattern = 'SDSS-I.fits'
    with urllib.request.urlopen('http://das.sdss.org/spectro/1d_26/0273/1d/spSpec-51957-0273-016.fit') as response:
        with tempfile.NamedTemporaryFile(prefix=sp_pattern) as tmp_file:
            shutil.copyfileobj(response, tmp_file)

            with warnings.catch_warnings():
                warnings.simplefilter('ignore', FITSFixedWarning)
                os.system(f'{compress} {tmp_file.name}')
                spec = Spectrum1D.read(tmp_file.name + ext[compress])

            assert isinstance(spec, Spectrum1D)
            assert spec.flux.size > 0
            assert spec.uncertainty.array.min() >= 0.0

            # Try again without compression suffix:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', FITSFixedWarning)
                os.system(f'mv {tmp_file.name}{ext[compress]} {tmp_file.name}')
                spec = Spectrum1D.read(tmp_file.name)

            assert isinstance(spec, Spectrum1D)
            assert spec.flux.size > 0
            assert spec.uncertainty.array.min() >= 0.0


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


@pytest.mark.parametrize("spectral_axis",
                         ['wavelength', 'frequency', 'energy', 'wavenumber'])
def test_tabular_fits_writer(tmpdir, spectral_axis):
    wlu = {'wavelength': u.AA, 'frequency': u.GHz, 'energy': u.eV,
           'wavenumber': u.cm**-1}
    # Create a small data set
    disp = np.arange(1,1.1,0.01)*wlu[spectral_axis]
    flux = np.ones(len(disp))*1.e-14*u.Jy
    unc = StdDevUncertainty(0.01*flux)
    spectrum = Spectrum1D(flux=flux, spectral_axis=disp, uncertainty=unc)
    tmpfile = str(tmpdir.join('_tst.fits'))
    spectrum.write(tmpfile, format='tabular-fits')

    # Read it in and check against the original
    spec = Spectrum1D.read(tmpfile)
    assert spec.flux.unit == spectrum.flux.unit
    assert spec.spectral_axis.unit == spectrum.spectral_axis.unit
    assert quantity_allclose(spec.spectral_axis, spectrum.spectral_axis)
    assert quantity_allclose(spec.flux, spectrum.flux)
    assert quantity_allclose(spec.uncertainty.quantity,
                             spectrum.uncertainty.quantity)

    # Test spectrum with different flux unit
    flux = np.random.normal(0., 1.e-9, disp.shape[0]) * u.W * u.m**-2 * u.AA**-1
    unc = StdDevUncertainty(0.1 * np.sqrt(np.abs(flux.value)) * flux.unit)
    spectrum = Spectrum1D(flux=flux, spectral_axis=disp, uncertainty=unc)

    # Try to overwrite the file
    with pytest.raises(OSError, match=r'File exists:'):
        spectrum.write(tmpfile, format='tabular-fits')
    spectrum.write(tmpfile, format='tabular-fits', overwrite=True)

    cmap = {spectral_axis: ('spectral_axis', wlu[spectral_axis]),
            'flux': ('flux', 'erg / (s cm**2 AA)'),
            'uncertainty': ('uncertainty', None)}

    # Read it back again and check against the original
    spec = Spectrum1D.read(tmpfile, format='tabular-fits', column_mapping=cmap)
    assert spec.flux.unit == u.Unit('erg / (s cm**2 AA)')
    assert spec.spectral_axis.unit == spectrum.spectral_axis.unit
    assert quantity_allclose(spec.spectral_axis, spectrum.spectral_axis)
    assert quantity_allclose(spec.flux, spectrum.flux)
    assert quantity_allclose(spec.uncertainty.quantity,
                             spectrum.uncertainty.quantity)


def test_tabular_fits_header(tmpdir):
    # Create a small data set + header with reserved FITS keywords
    disp = np.linspace(1, 1.2, 21) * u.AA
    flux = np.random.normal(0., 1.0e-14, disp.shape[0]) * u.Jy
    hdr = fits.header.Header({'TELESCOP': 'Leviathan', 'APERTURE': 1.8,
                              'OBSERVER': 'Parsons', 'NAXIS': 1, 'NAXIS1': 8})

    spectrum = Spectrum1D(flux=flux, spectral_axis=disp, meta={'header': hdr})
    tmpfile = str(tmpdir.join('_tst.fits'))
    spectrum.write(tmpfile, format='tabular-fits')

    # Read it in and check against the original
    hdulist = fits.open(tmpfile)
    assert hdulist[0].header['NAXIS'] == 0
    assert hdulist[1].header['NAXIS'] == 2
    assert hdulist[1].header['NAXIS2'] == disp.shape[0]
    assert hdulist[1].header['OBSERVER'] == 'Parsons'
    hdulist.close()

    # Now write with updated header information from spectrum.meta
    spectrum.meta.update({'OBSERVER': 'Rosse', 'EXPTIME': 32.1, 'NAXIS2': 12})
    spectrum.write(tmpfile, format='tabular-fits', overwrite=True,
                   update_header=True)

    hdulist = fits.open(tmpfile)
    assert hdulist[1].header['NAXIS2'] == disp.shape[0]
    assert hdulist[1].header['OBSERVER'] == 'Rosse'
    assert_allclose(hdulist[1].header['EXPTIME'], 3.21e1)
    hdulist.close()

    # Test that unsupported types (dict) are not added to written header
    spectrum.meta['MYHEADER'] = {'OBSDATE': '1848-02-26', 'TARGET': 'M51'}
    spectrum.write(tmpfile, format='tabular-fits', overwrite=True,
                   update_header=True)

    hdulist = fits.open(tmpfile)
    assert 'MYHEADER' not in hdulist[0].header
    assert 'MYHEADER' not in hdulist[1].header
    assert 'OBSDATE' not in hdulist[0].header
    assert 'OBSDATE' not in hdulist[1].header
    hdulist.close()


@pytest.mark.remote_data
def test_apstar_loader():
    '''Test remote read and automatic recognition of apStar spec from URL.
    '''
    apstar_url = ("https://data.sdss.org/sas/dr16/apogee/spectro/redux/r12/"
                  "stars/apo25m/N7789/apStar-r12-2M00005414+5522241.fits")
    spec = Spectrum1D.read(apstar_url)

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0
    assert spec.uncertainty.array.min() >= 0.0


@pytest.mark.remote_data
def test_apvisit_loader():
    '''Test remote read and automatic recognition of apvisit spec from URL.
    '''
    apvisit_url = ("https://data.sdss.org/sas/dr16/apogee/spectro/redux/r12/"
                   "visit/apo25m/N7789/5094/55874/"
                   "apVisit-r12-5094-55874-123.fits")
    spec = Spectrum1D.read(apvisit_url)

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0
    assert spec.uncertainty.array.min() >= 0.0


@pytest.mark.remote_data
def test_aspcapstar_loader():
    '''Test remote read and automatic recognition of aspcapStar spec from URL.
    '''
    aspcap_url = ("https://data.sdss.org/sas/dr16/apogee/spectro/aspcap/r12/"
                  "l33/apo25m/N7789/aspcapStar-r12-2M00005414+5522241.fits")
    spec = Spectrum1D.read(aspcap_url)

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0
    assert spec.uncertainty.array.min() >= 0.0


@pytest.mark.remote_data
def test_muscles_loader():
    '''Test remote read and automatic recognition of muscles spec from URL.
    '''
    url = ("https://archive.stsci.edu/missions/hlsp/muscles/gj1214/"
           "hlsp_muscles_multi_multi_gj1214_broadband_v22_const-res-sed.fits")
    spec = Spectrum1D.read(url)

    assert isinstance(spec, Spectrum1D)
    assert len(spec.flux) == len(spec.spectral_axis) > 50000
    assert spec.uncertainty.array.min() >= 0.0
