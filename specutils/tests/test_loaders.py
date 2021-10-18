import os
import sys
import shutil
import tempfile
import urllib
import warnings
import uuid
import pytest

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.table import Table
from astropy.units import UnitsWarning
from astropy.wcs import FITSFixedWarning, WCS
from astropy.io.registry import IORegistryError
from astropy.modeling import models
from astropy.tests.helper import quantity_allclose
from astropy.nddata import StdDevUncertainty

from numpy.testing import assert_allclose

from .conftest import remote_access
from .. import Spectrum1D, SpectrumCollection, SpectrumList
from ..io import get_loaders_by_extension
from ..io.default_loaders import subaru_pfs_spec
from ..io.default_loaders.sdss import _sdss_wcs_to_log_wcs

# NOTE: Python can be built without bz2 or lzma.
try:
    import bz2  # noqa
except ImportError:
    HAS_BZ2 = False
else:
    HAS_BZ2 = True

try:
    import lzma  # noqa
except ImportError:
    HAS_LZMA = False
else:
    HAS_LZMA = True


EBOSS_SPECTRUM_URL = 'https://data.sdss.org/sas/dr16/eboss/spectro/redux/v5_13_0/spectra/lite/4055/spec-4055-55359-0596.fits'


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
    wave = np.arange(1, 1.1, 0.01)*u.AA
    flux = np.ones(len(wave))*1.e-14*u.Jy
    uncertainty = 0.01*flux
    table = Table([wave, flux, uncertainty], names=["wave", "flux", "uncertainty"])
    tmpfile = str(tmpdir.join('_tst.ecsv'))
    table.write(tmpfile, format='ascii.ecsv')

    # Read it in and check against the original
    spectrum = Spectrum1D.read(tmpfile, format='ECSV')
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

    # HDUList case
    hdulist = fits.open(remote_data_path)
    spec = Spectrum1D.read(hdulist, format="HST/COS")
    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0


@remote_access([{'id': '1481192', 'filename': 'STIS_FUV.fits'},
                {'id': '1481185', 'filename': 'STIS_NUV.fits'},
                {'id': '1481183', 'filename': 'STIS_CCD.fits'}])
def test_hst_stis(remote_data_path):
    spec = Spectrum1D.read(remote_data_path)

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0

    # HDUList case
    hdulist = fits.open(remote_data_path)
    spec = Spectrum1D.read(hdulist, format="HST/STIS")
    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0


@pytest.mark.remote_data
def test_manga_cube():
    url = 'https://dr15.sdss.org/sas/dr15/manga/spectro/redux/v2_4_3/8485/stack/manga-8485-1901-LOGCUBE.fits.gz'
    spec = Spectrum1D.read(url, format='MaNGA cube')

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0
    assert spec.meta['header']['INSTRUME'] == 'MaNGA'
    assert spec.shape == (34, 34, 4563)


@pytest.mark.remote_data
def test_manga_rss():
    url = 'https://dr15.sdss.org/sas/dr15/manga/spectro/redux/v2_4_3/8485/stack/manga-8485-1901-LOGRSS.fits.gz'
    spec = Spectrum1D.read(url, format='MaNGA rss')

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0
    assert spec.meta['header']['INSTRUME'] == 'MaNGA'
    assert spec.shape == (171, 4563)


@pytest.mark.remote_data
def test_sdss_spec():
    sp_pattern = 'spec-4055-55359-0596.fits'
    with urllib.request.urlopen(EBOSS_SPECTRUM_URL) as response:
        # Read from open file object
        spec = Spectrum1D.read(response, format="SDSS-III/IV spec")
        assert isinstance(spec, Spectrum1D)
        assert spec.flux.size > 0

    with urllib.request.urlopen(EBOSS_SPECTRUM_URL) as response:
        # On Windows, NamedTemporaryFile cannot be opened a second time while
        #  already being open, so we avoid using that method.
        with tempfile.TemporaryDirectory() as tmp_dir:
            file_path = os.path.join(tmp_dir, sp_pattern)

            with open(file_path, 'wb') as tmp_file:
                shutil.copyfileobj(response, tmp_file)

                # Read from local disk via filename
                print(tmp_file.name)
                spec = Spectrum1D.read(tmp_file.name)

                assert isinstance(spec, Spectrum1D)
                assert spec.flux.size > 0

                # Read from HDUList object
                hdulist = fits.open(tmp_file.name)
                spec = Spectrum1D.read(hdulist)
                assert isinstance(spec, Spectrum1D)
                assert spec.flux.size > 0

                # Read from file handle
                fileio = open(tmp_file.name, mode='rb')
                spec = Spectrum1D.read(fileio)
                assert isinstance(spec, Spectrum1D)
                assert spec.flux.size > 0


@pytest.mark.remote_data
def test_sdss_spspec():
    sp_pattern = 'spSpec-51957-0273-016.fit'
    with urllib.request.urlopen('http://das.sdss.org/spectro/1d_26/0273/1d/spSpec-51957-0273-016.fit') as response:
        # Read from open file object
        spec = Spectrum1D.read(response, format="SDSS-I/II spSpec")
        assert isinstance(spec, Spectrum1D)
        assert spec.flux.size > 0

    with urllib.request.urlopen('http://das.sdss.org/spectro/1d_26/0273/1d/spSpec-51957-0273-016.fit') as response:
        # On Windows, NamedTemporaryFile cannot be opened a second time while
        #  already being open, so we avoid using that method.
        with tempfile.TemporaryDirectory() as tmp_dir:
            file_path = os.path.join(tmp_dir, sp_pattern)

            with open(file_path, 'wb') as tmp_file:
                shutil.copyfileobj(response, tmp_file)

                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', FITSFixedWarning)
                    spec = Spectrum1D.read(tmp_file.name)

                assert isinstance(spec, Spectrum1D)
                assert spec.flux.size > 0

                # Read from HDUList object
                hdulist = fits.open(tmp_file.name)
                spec = Spectrum1D.read(hdulist)
                assert isinstance(spec, Spectrum1D)
                assert spec.flux.size > 0

                # Read from file handle
                fileio = open(tmp_file.name, mode='rb')
                spec = Spectrum1D.read(fileio)
                assert isinstance(spec, Spectrum1D)
                assert spec.flux.size > 0


@pytest.mark.remote_data
def test_sdss_spec_stream():
    """Test direct read and recognition of SDSS-III/IV spec from remote URL,
    i.e. do not rely on filename pattern.
    """
    spec = Spectrum1D.read(EBOSS_SPECTRUM_URL)

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0
    assert spec.uncertainty.array.min() >= 0.0


@pytest.mark.remote_data
def test_sdss_spspec_stream():
    """Test direct read and recognition of SDSS-I/II spSpec from remote URL,
    i.e. do not rely on filename pattern.
    """
    sdss_url = 'http://das.sdss.org/spectro/1d_26/0273/1d/spSpec-51957-0273-016.fit'
    spec = Spectrum1D.read(sdss_url)

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0
    assert spec.uncertainty.array.min() >= 0.0


@pytest.mark.skipif('sys.platform.startswith("win")',
                    reason='Uncertain availability of compression utilities')
@pytest.mark.remote_data
@pytest.mark.parametrize('compress', ['gzip', 'bzip2', 'xz'])
def test_sdss_compressed(compress, tmp_path):
    """Test automatic recognition of supported compression formats.
    """
    ext = {'gzip': '.gz', 'bzip2': '.bz2', 'xz': '.xz'}
    if compress == 'bzip2' and not HAS_BZ2:
        pytest.xfail("Python installation has no bzip2 support")
    if compress == 'xz' and not HAS_LZMA:
        pytest.xfail("Python installation has no lzma support")

    # Deliberately not using standard filename pattern to test header info.
    tmp_filename = tmp_path / 'SDSS-I.fits'
    with urllib.request.urlopen('http://das.sdss.org/spectro/1d_26/0273/1d/spSpec-51957-0273-016.fit') as response:
        with open(tmp_filename, 'wb') as tmp_file:
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


@pytest.mark.remote_data
def test_sdss_spplate():
    """Test loading of multi-object spectrum from SDSS `spPlate` format FITS file.
    """
    with urllib.request.urlopen('http://das.sdss.org/spectro/1d_26/0273/1d/spSpec-51957-0273-016.fit') as response:
        # Read reference spectrum from open file object
        spec = Spectrum1D.read(response, format="SDSS-I/II spSpec")
        assert isinstance(spec, Spectrum1D)
        assert spec.flux.size > 0
        specid = spec.meta['header']['FIBERID']

    with urllib.request.urlopen('https://data.sdss.org/sas/dr8/sdss/spectro/redux/26/0273/spPlate-0273-51957.fits') as response:
        # Read "plate" spectrum with 2D flux array from open file object
        plate = Spectrum1D.read(response, format="SDSS spPlate")
        assert isinstance(plate, Spectrum1D)
        assert plate.flux.ndim == 2
        assert plate.flux.shape[0] == 640
        assert quantity_allclose(spec.spectral_axis, plate.spectral_axis)
        assert quantity_allclose(spec.flux, plate.flux[specid-1])

    with urllib.request.urlopen('https://data.sdss.org/sas/dr8/sdss/spectro/redux/26/0273/spPlate-0273-51957.fits') as response:
        with tempfile.NamedTemporaryFile() as tmp_file:
            shutil.copyfileobj(response, tmp_file)

            # Read from local disk via file signature
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', FITSFixedWarning)
                plate = Spectrum1D.read(tmp_file.name, limit=32)

            assert isinstance(plate, Spectrum1D)
            assert plate.flux.ndim == 2
            assert plate.flux.shape[0] == 32
            assert quantity_allclose(spec.spectral_axis, plate.spectral_axis)
            assert quantity_allclose(spec.flux, plate.flux[specid-1])

            # Read from HDUList object
            hdulist = fits.open(tmp_file.name)
            plate = Spectrum1D.read(hdulist, limit=32)
            assert plate.flux.shape[0] == 32
            assert quantity_allclose(spec.spectral_axis, plate.spectral_axis)
            assert quantity_allclose(spec.flux, plate.flux[specid-1])

            # Read from file handle
            fileio = open(tmp_file.name, mode='rb')
            plate = Spectrum1D.read(fileio, limit=32)
            assert plate.flux.shape[0] == 32
            assert quantity_allclose(spec.spectral_axis, plate.spectral_axis)
            assert quantity_allclose(spec.flux, plate.flux[specid-1])


@pytest.mark.parametrize("name", ['file.fit', 'file.fits', 'file.dat'])
def test_no_reader_matches(name):
    """If no reader matches a file, check that the correct error is raised.
    This test serves a second purpose: A badly written identifier
    function might raise an error as supposed to returning False when
    it cannot identify a file.  The fact that this test passes means
    that at the very least all identifier functions that have been
    tried for that file ending did not fail with an error.
    """
    with tempfile.TemporaryDirectory() as tmpdirname:
        filename = os.path.join(tmpdirname, name)
        with open(filename, 'w') as fp:
            fp.write('asdfadasdadvzxcv')

        with pytest.raises(IORegistryError):
            spec = Spectrum1D.read(filename)


@remote_access([{'id': '3359174', 'filename': 'linear_fits_solution.fits'}])
def test_iraf_linear(remote_data_path):

    spectrum_1d = Spectrum1D.read(remote_data_path, format='iraf')

    assert isinstance(spectrum_1d, Spectrum1D)
    assert quantity_allclose(spectrum_1d.wavelength[0],
                             u.Quantity(3514.56625402, unit='Angstrom'))
    assert quantity_allclose(spectrum_1d.wavelength[100],
                             u.Quantity(3514.56625402, unit='Angstrom') +
                             u.Quantity(0.653432383823 * 100, unit='Angstrom'))


@remote_access([{'id': '3359180', 'filename': 'log-linear_fits_solution.fits'}])
def test_iraf_log_linear(remote_data_path):

    with pytest.raises(NotImplementedError):
        assert Spectrum1D.read(remote_data_path, format='iraf')


@remote_access([{'id': '3359190', 'filename': 'non-linear_fits_solution_cheb.fits'}])
def test_iraf_non_linear_chebyshev(remote_data_path):
    chebyshev_model = models.Chebyshev1D(degree=2, domain=[1616, 3259])
    chebyshev_model.c0.value = 5115.64008186
    chebyshev_model.c1.value = 535.515983712
    chebyshev_model.c2.value = -0.779265625182

    wavelength_axis = chebyshev_model(range(1, 4097)) * u.angstrom

    spectrum_1d = Spectrum1D.read(remote_data_path, format='iraf')

    assert isinstance(spectrum_1d, Spectrum1D)
    assert_allclose(wavelength_axis, spectrum_1d.wavelength)

    # Read from HDUList
    hdulist = fits.open(remote_data_path)
    spectrum_1d = Spectrum1D.read(hdulist, format='iraf')
    assert isinstance(spectrum_1d, Spectrum1D)
    assert_allclose(wavelength_axis, spectrum_1d.wavelength)


@remote_access([{'id': '3359194', 'filename': 'non-linear_fits_solution_legendre.fits'}])
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


@remote_access([{'id': '3359196', 'filename': 'non-linear_fits_solution_linear-spline.fits'}])
def test_iraf_non_linear_linear_spline(remote_data_path):

    with pytest.raises(NotImplementedError):
        assert Spectrum1D.read(remote_data_path, format='iraf')


@remote_access([{'id': '3359200', 'filename': 'non-linear_fits_solution_cubic-spline.fits'}])
def test_iraf_non_linear_cubic_spline(remote_data_path):

    with pytest.raises(NotImplementedError):
        assert Spectrum1D.read(remote_data_path, format='iraf')


@pytest.mark.remote_data
def test_iraf_multispec_chebyshev():
    """Test loading of SpectrumCollection from IRAF MULTISPEC format FITS file -
    nonlinear 2D WCS with Chebyshev solution (FTYPE=1).
    """
    iraf_url = 'https://github.com/astropy/specutils/raw/legacy-specutils/specutils/io/tests/files'
    # Read reference ASCII spectrum from remote file.
    spec10 = Table.read(iraf_url + '/AAO_11.txt', format='ascii.no_header',
                        data_start=175, names=['spectral_axis', 'flux'])

    # Read full collection of 51 spectra from open FITS file object
    speccol = SpectrumCollection.read(iraf_url + '/AAO.fits')

    assert len(speccol) == 51
    # Numpy allclose does support quantities, but not with unit 'adu'!
    assert quantity_allclose(speccol[10].spectral_axis, spec10['spectral_axis'] * u.AA)
    assert quantity_allclose(speccol[10].flux, spec10['flux'] * u.adu)


@pytest.mark.remote_data
def test_iraf_multispec_legendre():
    """Test loading of SpectrumCollection from IRAF MULTISPEC format FITS file -
    nonlinear 2D WCS with Legendre solution (FTYPE=2).
    """
    iraf_url = 'https://github.com/astropy/specutils/raw/legacy-specutils/specutils/io/tests/files'
    # Read reference ASCII spectrum from remote file.
    spec10 = Table.read(iraf_url + '/TRES.dat', format='ascii.no_header',
                        data_start=127, names=['spectral_axis', 'flux'])

    # Read full collection of 51 spectra from remote FITS file object
    speccol = SpectrumCollection.read(iraf_url + '/TRES.fits')

    assert len(speccol) == 51
    # The reference spectrum flux is normalised/flatfielded, just check the wavelength solution
    assert_allclose(speccol[10].spectral_axis, spec10['spectral_axis'] * u.AA)


@pytest.mark.parametrize("spectral_axis",
                         ['wavelength', 'frequency', 'energy', 'wavenumber'])
def test_tabular_fits_writer(tmpdir, spectral_axis):
    wlu = {'wavelength': u.AA, 'frequency': u.GHz, 'energy': u.eV,
           'wavenumber': u.cm**-1}
    # Create a small data set
    disp = np.arange(1, 1.1, 0.01) * wlu[spectral_axis]
    flux = np.ones(len(disp)) * 1.e-14 * u.Jy
    unc = StdDevUncertainty(0.01 * flux)
    if spectral_axis not in ('wavelength', ):
        disp = np.flip(disp)

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
    with pytest.raises(OSError, match=r'File .*exists'):
        spectrum.write(tmpfile, format='tabular-fits')
    spectrum.write(tmpfile, format='tabular-fits', overwrite=True)

    # Map to alternative set of units
    cmap = {spectral_axis: ('spectral_axis', 'micron'),
            'flux': ('flux', 'erg / (s cm**2 AA)'),
            'uncertainty': ('uncertainty', None)}

    # Read it back again and check against the original
    spec = Spectrum1D.read(tmpfile, format='tabular-fits', column_mapping=cmap)
    assert spec.flux.unit == u.Unit('erg / (s cm**2 AA)')
    assert spec.spectral_axis.unit == u.um
    assert quantity_allclose(spec.spectral_axis, spectrum.spectral_axis)
    assert quantity_allclose(spec.flux, spectrum.flux)
    assert quantity_allclose(spec.uncertainty.quantity,
                             spectrum.uncertainty.quantity)


@pytest.mark.parametrize("ndim", range(1, 4))
@pytest.mark.parametrize("spectral_axis",
                         ['wavelength', 'frequency', 'energy', 'wavenumber'])
def test_tabular_fits_multid(tmpdir, ndim, spectral_axis):
    wlu = {'wavelength': u.AA, 'frequency': u.GHz, 'energy': u.eV,
           'wavenumber': u.cm**-1}
    # Create a small data set with ndim-D flux + uncertainty
    disp = np.arange(1, 1.1, 0.01) * wlu[spectral_axis]
    shape = (3, 2, 4)[:ndim+1] + disp.shape
    flux = np.random.normal(0., 1.e-9, shape) * u.W * u.m**-2 * u.AA**-1
    unc = StdDevUncertainty(0.01 * np.random.sample(shape))
    if spectral_axis not in ('wavelength', ):
        disp = np.flip(disp)

    spectrum = Spectrum1D(flux=flux, spectral_axis=disp, uncertainty=unc)
    tmpfile = str(tmpdir.join('_tst.fits'))
    spectrum.write(tmpfile, format='tabular-fits')

    # Read it in and check against the original
    spec = Spectrum1D.read(tmpfile)
    assert spec.flux.unit == spectrum.flux.unit
    assert spec.spectral_axis.unit == spectrum.spectral_axis.unit
    assert spec.flux.shape == flux.shape
    assert spec.uncertainty.array.shape == flux.shape
    assert quantity_allclose(spec.spectral_axis, spectrum.spectral_axis)
    assert quantity_allclose(spec.flux, spectrum.flux)
    assert quantity_allclose(spec.uncertainty.quantity,
                             spectrum.uncertainty.quantity)

    # Test again, using `column_mapping` to convert to different spectral axis and flux units
    cmap = {spectral_axis: ('spectral_axis', 'THz'),
            'flux': ('flux', 'erg / (s cm**2 AA)'),
            'uncertainty': ('uncertainty', None)}

    spec = Spectrum1D.read(tmpfile, format='tabular-fits', column_mapping=cmap)
    assert spec.flux.unit == u.Unit('erg / (s cm**2 AA)')
    assert spec.spectral_axis.unit == u.THz
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


def test_tabular_fits_autowrite(tmpdir):
    """Test writing of Spectrum1D with automatic selection of BINTABLE format."""
    disp = np.linspace(1, 1.2, 21) * u.AA
    flux = np.random.normal(0., 1.0e-14, disp.shape[0]) * u.W / (u.m**2 * u.AA)
    hdr = fits.header.Header({'TELESCOP': 'Leviathan', 'APERTURE': 1.8,
                              'OBSERVER': 'Parsons', 'NAXIS': 1, 'NAXIS1': 8})

    spectrum = Spectrum1D(flux=flux, spectral_axis=disp, meta={'header': hdr})
    tmpfile = str(tmpdir.join('_tst.fits'))
    spectrum.write(tmpfile)

    # Read it in and check against the original
    with fits.open(tmpfile) as hdulist:
        assert hdulist[0].header['NAXIS'] == 0
        assert hdulist[1].header['NAXIS'] == 2
        assert hdulist[1].header['NAXIS2'] == disp.shape[0]

    # Trigger exception for illegal HDU (primary HDU only accepts IMAGE_HDU)
    with pytest.raises(ValueError, match=r'FITS does not support BINTABLE'):
        spectrum.write(tmpfile, format='tabular-fits', overwrite=True, hdu=0)

    # Test automatic selection of wcs1d format, which will fail without suitable wcs
    with pytest.raises(ValueError, match=r'Only Spectrum1D objects with valid WCS'):
        spectrum.write(tmpfile, overwrite=True, hdu=0)

    tmpfile = str(tmpdir.join('_wcs.fits'))
    with pytest.raises(ValueError, match=r'Only Spectrum1D objects with valid WCS'):
        spectrum.write(tmpfile, overwrite=True)


@pytest.mark.skipif('sys.platform.startswith("win")',
                    reason='Uncertain availability of compression utilities')
@pytest.mark.parametrize('compress', ['gzip', 'bzip2', 'xz'])
def test_tabular_fits_compressed(compress, tmpdir):
    """Test automatic recognition of supported compression formats for BINTABLE.
    """
    ext = {'gzip': '.gz', 'bzip2': '.bz2', 'xz': '.xz'}
    if compress == 'bzip2' and not HAS_BZ2:
        pytest.xfail("Python installation has no bzip2 support")
    if compress == 'xz' and not HAS_LZMA:
        pytest.xfail("Python installation has no lzma support")

    # Create a small data set
    disp = np.linspace(1, 1.2, 23) * u.AA
    flux = np.random.normal(0., 1.0e-14, disp.shape[0]) * u.Jy
    unc = StdDevUncertainty(0.01 * np.abs(flux))

    spectrum = Spectrum1D(flux=flux, spectral_axis=disp, uncertainty=unc)
    tmpfile = str(tmpdir.join('_tst.fits'))
    spectrum.write(tmpfile, format='tabular-fits')

    # Deliberately not using standard filename pattern to test header info.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FITSFixedWarning)
        os.system(f'{compress} {tmpfile}')
        spec = Spectrum1D.read(tmpfile + ext[compress])

    assert isinstance(spec, Spectrum1D)
    assert spec.spectral_axis.shape[0] == len(disp)
    assert spec.flux.size == len(disp)
    assert spec.uncertainty.array.min() >= 0.0
    assert quantity_allclose(spec.flux, spectrum.flux)

    # Try again without compression suffix:
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FITSFixedWarning)
        os.system(f'mv {tmpfile}{ext[compress]} {tmpfile}')
        spec = Spectrum1D.read(tmpfile)

    assert isinstance(spec, Spectrum1D)
    assert spec.spectral_axis.shape[0] == len(disp)
    assert spec.flux.size == len(disp)
    assert spec.uncertainty.array.min() >= 0.0
    assert quantity_allclose(spec.flux, spectrum.flux)


@pytest.mark.parametrize("spectral_axis",
                         ['wavelength', 'frequency', 'energy', 'wavenumber'])
def test_wcs1d_fits_writer(tmpdir, spectral_axis):
    """Test write/read for Spectrum1D with WCS-constructed spectral_axis."""
    wlunits = {'wavelength': 'Angstrom', 'frequency': 'GHz', 'energy': 'eV',
               'wavenumber': 'cm**-1'}
    # Header dictionary for constructing WCS
    hdr = {'CTYPE1': spectral_axis, 'CUNIT1': wlunits[spectral_axis],
           'CRPIX1': 1, 'CRVAL1': 1, 'CDELT1': 0.01}
    # Create a small data set
    flux = np.arange(1, 11)**2 * 1.e-14 * u.Jy
    wlu = u.Unit(hdr['CUNIT1'])
    wl0 = hdr['CRVAL1']
    dwl = hdr['CDELT1']
    disp = np.arange(wl0, wl0 + (len(flux) - 0.5) * dwl, dwl) * wlu

    spectrum = Spectrum1D(flux=flux, wcs=WCS(hdr))
    tmpfile = str(tmpdir.join('_tst.fits'))
    spectrum.write(tmpfile, hdu=0)

    # Read it in and check against the original
    spec = Spectrum1D.read(tmpfile)
    assert spec.flux.unit == spectrum.flux.unit
    assert spec.spectral_axis.unit == spectrum.spectral_axis.unit
    assert quantity_allclose(spec.spectral_axis, spectrum.spectral_axis)
    assert quantity_allclose(spec.spectral_axis, disp)
    assert quantity_allclose(spec.flux, spectrum.flux)

    # Read from HDUList
    hdulist = fits.open(tmpfile)
    spec = Spectrum1D.read(hdulist, format='wcs1d-fits')
    assert isinstance(spec, Spectrum1D)
    assert quantity_allclose(spec.spectral_axis, spectrum.spectral_axis)
    assert quantity_allclose(spec.flux, spectrum.flux)


@pytest.mark.parametrize("hdu", range(3))
def test_wcs1d_fits_hdus(tmpdir, hdu):
    """Test writing of Spectrum1D in WCS1D format to different IMAGE_HDUs."""
    # Header dictionary for constructing WCS
    hdr = {'CTYPE1': 'wavelength', 'CUNIT1': 'um',
           'CRPIX1': 1, 'CRVAL1': 1, 'CDELT1': 0.01}
    # Create a small data set
    flu = u.W / (u.m**2 * u.nm)
    flux = np.arange(1, 11)**2 * 1.e-14 * flu

    spectrum = Spectrum1D(flux=flux, wcs=WCS(hdr))
    tmpfile = str(tmpdir.join('_tst.fits'))
    spectrum.write(tmpfile, hdu=hdu, format='wcs1d-fits')

    # Read it in and check against the original
    with fits.open(tmpfile) as hdulist:
        assert hdulist[hdu].is_image
        assert hdulist[hdu].header['NAXIS'] == 1
        assert hdulist[hdu].header['NAXIS1'] == flux.shape[0]
        assert u.Unit(hdulist[hdu].header['CUNIT1']) == u.Unit(hdr['CUNIT1'])
        assert quantity_allclose(hdulist[hdu].data * flu, flux)

    # Test again with automatic format selection by filename pattern
    tmpfile = str(tmpdir.join('_wcs.fits'))
    spectrum.write(tmpfile, hdu=hdu)
    with fits.open(tmpfile) as hdulist:
        assert hdulist[hdu].is_image
        assert quantity_allclose(hdulist[hdu].data * flu, flux)


@pytest.mark.parametrize("spectral_axis",
                         ['wavelength', 'frequency', 'energy', 'wavenumber'])
def test_wcs1d_fits_multid(tmpdir, spectral_axis):
    """Test spectrum with WCS-1D spectral_axis and higher dimension in flux."""
    wlunits = {'wavelength': 'Angstrom', 'frequency': 'GHz', 'energy': 'eV',
               'wavenumber': 'cm**-1'}
    # Header dictionary for constructing WCS
    hdr = {'CTYPE1': spectral_axis, 'CUNIT1': wlunits[spectral_axis],
           'CRPIX1': 1, 'CRVAL1': 1, 'CDELT1': 0.01}
    # Create a small data set
    flux = np.arange(1, 11)**2 * 1.e-14 * u.Jy
    wlu = u.Unit(hdr['CUNIT1'])
    wl0 = hdr['CRVAL1']
    dwl = hdr['CDELT1']
    disp = np.arange(wl0, wl0 + len(flux[1:]) * dwl, dwl) * wlu

    # Construct up to 4D flux array, write and read (no auto-identify)
    shape = [-1, 1]
    for i in range(2, 5):
        flux = flux * np.arange(i, i+5).reshape(*shape)
        spectrum = Spectrum1D(flux=flux, wcs=WCS(hdr))
        tmpfile = str(tmpdir.join(f'_{i}d.fits'))
        spectrum.write(tmpfile, format='wcs1d-fits')

        spec = Spectrum1D.read(tmpfile, format='wcs1d-fits')
        assert spec.flux.ndim == i
        assert quantity_allclose(spec.spectral_axis, disp)
        assert quantity_allclose(spec.spectral_axis, spectrum.spectral_axis)
        assert quantity_allclose(spec.flux, spectrum.flux)
        shape.append(1)

    # Test exception for NAXIS > 4
    flux = flux * np.arange(i+1, i+6).reshape(*shape)
    spectrum = Spectrum1D(flux=flux, wcs=WCS(hdr))
    tmpfile = str(tmpdir.join(f'_{i+1}d.fits'))
    spectrum.write(tmpfile, format='wcs1d-fits')

    with pytest.raises(ValueError, match='input to wcs1d_fits_loader is > 4D'):
        spec = Spectrum1D.read(tmpfile, format='wcs1d-fits')


@pytest.mark.parametrize("spectral_axis", ['wavelength', 'frequency'])
def test_wcs1d_fits_non1d(tmpdir, spectral_axis):
    """Test exception on trying to load FITS with 2D flux and irreducible WCS
    spectral_axis.
    """
    wlunits = {'wavelength': 'Angstrom', 'frequency': 'GHz', 'energy': 'eV',
               'wavenumber': 'cm**-1'}
    # Header dictionary for constructing WCS
    hdr = {'CTYPE1': spectral_axis, 'CUNIT1': wlunits[spectral_axis],
           'CRPIX1': 1, 'CRVAL1': 1, 'CDELT1': 0.01}
    # Create a small 2D data set
    flux = np.arange(1, 11)**2 * np.arange(4).reshape(-1, 1) * 1.e-14 * u.Jy
    spectrum = Spectrum1D(flux=flux, wcs=WCS(hdr))
    tmpfile = str(tmpdir.join(f'_{2}d.fits'))
    spectrum.write(tmpfile, format='wcs1d-fits')

    # Reopen file and update header with off-diagonal element
    hdulist = fits.open(tmpfile, mode='update')
    hdulist[0].header.update([('PC1_2', 0.2)])
    hdulist.close()

    with pytest.raises(ValueError, match='WCS cannot be reduced to 1D'):
        Spectrum1D.read(tmpfile, format='wcs1d-fits')


@pytest.mark.skipif('sys.platform.startswith("win")',
                    reason='Uncertain availability of compression utilities')
@pytest.mark.parametrize('compress', ['gzip', 'bzip2', 'xz'])
def test_wcs1d_fits_compressed(compress, tmpdir):
    """Test automatic recognition of supported compression formats for IMAGE/WCS.
    """
    ext = {'gzip': '.gz', 'bzip2': '.bz2', 'xz': '.xz'}
    if compress == 'bzip2' and not HAS_BZ2:
        pytest.xfail("Python installation has no bzip2 support")
    if compress == 'xz' and not HAS_LZMA:
        pytest.xfail("Python installation has no lzma support")

    # Header dictionary for constructing WCS
    hdr = {'CTYPE1': 'wavelength', 'CUNIT1': 'Angstrom',
           'CRPIX1': 1, 'CRVAL1': 1, 'CDELT1': 0.01}
    # Create a small data set
    flux = np.arange(1, 43)**2 * 1.e-14 * u.Jy
    wlu = u.Unit(hdr['CUNIT1'])
    wl0 = hdr['CRVAL1']
    dwl = hdr['CDELT1']
    disp = np.arange(wl0, wl0 + (len(flux) - 0.5) * dwl, dwl) * wlu

    spectrum = Spectrum1D(flux=flux, wcs=WCS(hdr))
    tmpfile = str(tmpdir.join('_tst.fits'))
    spectrum.write(tmpfile, hdu=0)

    # Deliberately not using standard filename pattern to test header info.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FITSFixedWarning)
        os.system(f'{compress} {tmpfile}')
        spec = Spectrum1D.read(tmpfile + ext[compress])

    assert isinstance(spec, Spectrum1D)
    assert quantity_allclose(spec.spectral_axis, disp)
    assert quantity_allclose(spec.flux, spectrum.flux)

    # Try again without compression suffix:
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FITSFixedWarning)
        os.system(f'mv {tmpfile}{ext[compress]} {tmpfile}')
        spec = Spectrum1D.read(tmpfile)

    assert isinstance(spec, Spectrum1D)
    assert quantity_allclose(spec.spectral_axis, disp)
    assert quantity_allclose(spec.flux, spectrum.flux)


@pytest.mark.remote_data
def test_apstar_loader():
    """Test remote read and automatic recognition of apStar spec from URL.
    """
    apstar_url = ("https://data.sdss.org/sas/dr16/apogee/spectro/redux/r12/"
                  "stars/apo25m/N7789/apStar-r12-2M00005414+5522241.fits")
    spec = Spectrum1D.read(apstar_url)

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0
    assert spec.flux.unit == 1e-17 * u.erg / (u.s * u.cm**2 * u.AA)
    assert spec.uncertainty.array.min() >= 0.0


@pytest.mark.remote_data
def test_apvisit_loader():
    """Test remote read and automatic recognition of apvisit spec from URL.
    """
    apvisit_url = ("https://data.sdss.org/sas/dr16/apogee/spectro/redux/r12/"
                   "visit/apo25m/N7789/5094/55874/"
                   "apVisit-r12-5094-55874-123.fits")
    spec = Spectrum1D.read(apvisit_url)

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0
    assert spec.flux.unit == 1e-17 * u.erg / (u.s * u.cm**2 * u.AA)
    assert spec.uncertainty.array.min() >= 0.0


@pytest.mark.remote_data
def test_aspcapstar_loader():
    """Test remote read and automatic recognition of aspcapStar spec from URL.
    """
    aspcap_url = ("https://data.sdss.org/sas/dr16/apogee/spectro/aspcap/r12/"
                  "l33/apo25m/N7789/aspcapStar-r12-2M00005414+5522241.fits")
    spec = Spectrum1D.read(aspcap_url)

    assert isinstance(spec, Spectrum1D)
    assert spec.flux.size > 0
    assert spec.uncertainty.array.min() >= 0.0


@pytest.mark.remote_data
def test_muscles_loader():
    """Test remote read and automatic recognition of muscles spec from URL.
    """
    url = ("https://archive.stsci.edu/missions/hlsp/muscles/gj1214/"
           "hlsp_muscles_multi_multi_gj1214_broadband_v22_const-res-sed.fits")
    spec = Spectrum1D.read(url)

    assert isinstance(spec, Spectrum1D)
    assert len(spec.flux) == len(spec.spectral_axis) > 50000
    assert spec.uncertainty.array.min() >= 0.0
    assert spec.spectral_axis.unit == u.AA
    assert spec.flux.unit == u.erg / (u.s * u.cm**2 * u.AA)

    # Read HDUList
    hdulist = fits.open(url)
    spec = Spectrum1D.read(hdulist, format="MUSCLES SED")
    assert isinstance(spec, Spectrum1D)


@pytest.mark.remote_data
def test_subaru_pfs_loader(tmpdir):
    """Test remote read and automatic recognition of Subaru PFS spec from URL.
    """
    pfs = "pfsObject-00000-0,0-000-00000001-01-0x395428ab.fits"
    url = f"https://github.com/Subaru-PFS/datamodel/raw/master/examples/{pfs}"

    assert subaru_pfs_spec.identify_pfs_spec(url, url)

    # PFS loader parses metadata from filename, cannot read directly from url
    tmpfile = str(tmpdir.join(pfs))
    with urllib.request.urlopen(url) as response:
        shutil.copyfileobj(response, open(tmpfile, mode='wb'))

    assert subaru_pfs_spec.identify_pfs_spec(pfs, open(tmpfile, mode='rb'))
    spec = Spectrum1D.read(tmpfile, format='Subaru-pfsObject')
    assert isinstance(spec, Spectrum1D)

    spec = Spectrum1D.read(tmpfile)
    assert isinstance(spec, Spectrum1D)
    assert len(spec.flux) == len(spec.spectral_axis) > 10000
    assert spec.spectral_axis.unit == u.nm
    assert spec.flux.unit == u.nJy


@remote_access([{'id': '3733958', 'filename': '1D-c0022498-344732.fits'}])
def test_spectrum1d_6dfgs_tabular(remote_data_path):
    spec = Spectrum1D.read(remote_data_path)

    assert spec.spectral_axis.unit == u.Unit("Angstrom")
    assert spec.flux.unit == u.Unit("count/s")

    # Read from HDUList object
    hdulist = fits.open(remote_data_path)
    spec = Spectrum1D.read(hdulist, format="6dFGS-tabular")
    assert isinstance(spec, Spectrum1D)
    assert spec.flux.unit == u.Unit("count/s")
    assert spec.flux.size > 0
    hdulist.close()


@remote_access([{'id': '3733958', 'filename': 'all-c0022498-344732v_spectrum0.fits'}])
def test_spectrum1d_6dfgs_split_v(remote_data_path):
    spec = Spectrum1D.read(remote_data_path)

    assert spec.spectral_axis.unit == u.Unit("Angstrom")
    assert spec.flux.unit == u.Unit("count/Angstrom")

    # Read from HDUList object
    hdulist = fits.open(remote_data_path)
    spec = Spectrum1D.read(hdulist, format="6dFGS-split")
    assert isinstance(spec, Spectrum1D)
    assert spec.flux.unit == u.Unit("count/Angstrom")
    assert spec.flux.size > 0
    hdulist.close()


@remote_access([{'id': '3733958', 'filename': 'all-c0022498-344732r_spectrum0.fits'}])
def test_spectrum1d_6dfgs_split_r(remote_data_path):
    spec = Spectrum1D.read(remote_data_path)

    assert spec.spectral_axis.unit == u.Unit("Angstrom")
    assert spec.flux.unit == u.Unit("count/Angstrom")

    # Read from HDUList object
    hdulist = fits.open(remote_data_path)
    spec = Spectrum1D.read(hdulist, format="6dFGS-split")
    assert isinstance(spec, Spectrum1D)
    assert spec.flux.unit == u.Unit("count/Angstrom")
    assert spec.flux.size > 0
    hdulist.close()


@remote_access([{'id': '3733958', 'filename': 'all-c0022498-344732combined_spectrum0.fits'}])
def test_spectrum1d_6dfgs_split_combined(remote_data_path):
    spec = Spectrum1D.read(remote_data_path)

    assert spec.spectral_axis.unit == u.Unit("Angstrom")
    assert spec.flux.unit == u.Unit("count/Angstrom")

    # Read from HDUList object
    hdulist = fits.open(remote_data_path)
    spec = Spectrum1D.read(hdulist, format="6dFGS-split")
    assert isinstance(spec, Spectrum1D)
    assert spec.flux.unit == u.Unit("count/Angstrom")
    assert spec.flux.size > 0
    hdulist.close()


@remote_access([{'id': '3733958', 'filename': 'all-c0022498-344732.fits'}])
def test_spectrum1d_6dfgs_combined(remote_data_path):
    specs = SpectrumList.read(remote_data_path)

    for spec in specs:
        assert spec.spectral_axis.unit == u.Unit("Angstrom")
        assert spec.flux.unit == u.Unit("count/Angstrom")

    assert len(specs) == 3

    # Read from HDUList object
    hdulist = fits.open(remote_data_path)
    specs = SpectrumList.read(hdulist, format="6dFGS-combined")
    for spec in specs:
        assert isinstance(spec, Spectrum1D)
        assert spec.flux.unit == u.Unit("count/Angstrom")
        assert spec.flux.size > 0
        assert spec.meta["sky"].flux.unit == u.Unit("count/Angstrom")
        assert spec.meta["sky"].flux.size > 0

    assert len(specs) == 3

    hdulist.close()


# Commented out until science only is discussed
# @pytest.mark.remote_data
# def test_2slaq_lrg_loader_science_only():
#     """Test remote read and automatic recognition of 2SLAQ-LRG data from URL.
#     """
#     url = ("https://datacentral.org.au/services/sov/81480/download/"
#            "gama.dr2.spectra.2slaq-lrg.spectrum_1d/J143529.78-004306.4_1.fit/")
#     spec = Spectrum1D.read(url)
#
#     assert spec.spectral_axis.unit == u.AA
#     assert spec.flux.unit == u.count / u.s
#     assert spec.uncertainty is None


@remote_access([{'id': '3970324', 'filename': 'J143529.78-004306.4_1.fit'}])
def test_2slaq_lrg_loader_science_and_sky(remote_data_path):
    """Test remote read and automatic recognition of 2SLAQ-LRG data from URL.
    """
    science, sky = SpectrumList.read(remote_data_path)

    assert science.spectral_axis.unit == u.AA
    assert science.flux.unit == u.count / u.s
    assert science.uncertainty is None

    assert sky.spectral_axis.unit == u.AA
    assert sky.flux.unit == u.count / u.s
    assert sky.uncertainty is None


@remote_access([
    {'id': '3895436', 'filename': '000002.fits'},
    {'id': '3895436', 'filename': '000003.fits'},
    {'id': '3895436', 'filename': '000004.fits'},
])
def test_spectrum_list_2dfgrs_single(remote_data_path):
    specs = SpectrumList.read(remote_data_path)

    assert len(specs) == 1

    for spec in specs:
        assert spec.spectral_axis.unit == u.Unit("Angstrom")


    # Read from HDUList object
    hdulist = fits.open(remote_data_path)
    specs = SpectrumList.read(hdulist, format="2dFGRS")
    for spec in specs:
        assert isinstance(spec, Spectrum1D)
        assert spec.spectral_axis.unit == u.Unit("Angstrom")

    assert len(specs) == 1

    hdulist.close()


@remote_access([{'id': '3895436', 'filename': '000001.fits'}])
def test_spectrum_list_2dfgrs_multiple(remote_data_path):
    specs = SpectrumList.read(remote_data_path)

    assert len(specs) == 2

    for spec in specs:
        assert spec.spectral_axis.unit == u.Unit("Angstrom")

    # Read from HDUList object
    hdulist = fits.open(remote_data_path)
    specs = SpectrumList.read(hdulist, format="2dFGRS")
    for spec in specs:
        assert isinstance(spec, Spectrum1D)
        assert spec.spectral_axis.unit == u.Unit("Angstrom")

    assert len(specs) == 2

    hdulist.close()


def test_sdss_wcs_handler():
    sdss_wcs = WCS(naxis=2)
    sdss_wcs.wcs.crval[0] = 3.57880000000000E+00
    sdss_wcs.wcs.cd = [[1.00000000000000E-04, 0], [0, 1]]
    sdss_wcs.wcs.cunit[0] = u.Unit('Angstrom')
    fixed_wcs = _sdss_wcs_to_log_wcs(sdss_wcs)
    dropped_sdss_wcs = sdss_wcs.dropaxis(1)
    dropped_sdss_wcs.wcs.cunit[0] = ''  # Cannot handle units in powers
    sdss_wave = 10 ** dropped_sdss_wcs.pixel_to_world(np.arange(10)) * u.Unit('Angstrom')
    fixed_wave = fixed_wcs.pixel_to_world(np.arange(10))
    assert quantity_allclose(sdss_wave, fixed_wave)


class TestAAOmega2dF:
    @remote_access([{'id': '4460981', 'filename': "OBJ0039red.fits"}])
    def test_with_rwss(self, remote_data_path):
        spectra = SpectrumList.read(
            remote_data_path, format="Data Central AAOmega",
        )
        assert len(spectra) == 139
        for spec in spectra:
            assert spec.meta.get("label") is not None
            assert spec.meta.get("header") is not None
            assert spec.meta.get("purpose") is not None
            assert spec.meta.get("fibre_index") is not None

    @remote_access([{'id': '4460981', 'filename': "OBJ0032red.fits"}])
    def test_without_rwss(self, remote_data_path):
        spectra = SpectrumList.read(
            remote_data_path, format="Data Central AAOmega",
        )
        assert len(spectra) == 153
        for spec in spectra:
            assert spec.meta.get("label") is not None
            assert spec.meta.get("header") is not None
            assert spec.meta.get("purpose") is not None
            assert spec.meta.get("fibre_index") is not None

    @remote_access([{'id': '4460981', 'filename': "OBJ0039red.fits"}])
    def test_with_rwss_guess(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path)
        assert len(spectra) == 139
        for spec in spectra:
            assert spec.meta.get("label") is not None
            assert spec.meta.get("header") is not None
            assert spec.meta.get("purpose") is not None
            assert spec.meta.get("fibre_index") is not None

    @remote_access([{'id': '4460981', 'filename': "OBJ0032red.fits"}])
    def test_without_rwss_guess(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path)
        assert len(spectra) == 153
        for spec in spectra:
            assert spec.meta.get("label") is not None
            assert spec.meta.get("header") is not None
            assert spec.meta.get("purpose") is not None
            assert spec.meta.get("fibre_index") is not None


@remote_access([
    {'id': "4460981", 'filename':"1812260046012353.fits"}, # 4 exts
    {'id': "4460981", 'filename':"1311160005010021.fits"} # 5 exts
])
def test_galah(remote_data_path):
    spectra = SpectrumList.read(remote_data_path, format="GALAH")
    # Should be main spectra, without sky, and normalised (not in 4 ext)
    nspec = len(spectra)
    if spectra[0].meta["galah_hdu_format"] == 4:
        assert nspec == 2
    elif spectra[0].meta["galah_hdu_format"] == 5:
        assert nspec == 3
    else:
        assert False, "Unknown format"
    # normalised
    if nspec == 3:
        assert spectra[0].flux.unit == u.Unit('') # dimensionless
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[0].uncertainty is None
        assert spectra[0].meta.get("label") == "normalised spectra"
        assert spectra[0].meta.get("header") is not None

        # drop the normalised spectra, so 4 and 5 should now look the same
        spectra = spectra[1:]

    # main spectra
    assert spectra[0].flux.unit == u.count
    assert spectra[0].spectral_axis.unit == u.Angstrom
    assert isinstance(spectra[0].uncertainty, StdDevUncertainty)
    assert spectra[0].meta.get("label") is not None
    assert spectra[0].meta.get("header") is not None

    # No sky
    assert spectra[1].spectral_axis.unit == u.Angstrom
    assert spectra[1].flux.unit == u.count
    assert isinstance(spectra[1].uncertainty, StdDevUncertainty)
    assert spectra[1].meta.get("label") is not None
    assert spectra[1].meta.get("header") is not None


@pytest.mark.xfail(reason="Format is ambiguous")
@remote_access([
    {'id': "4460981", 'filename':"1812260046012353.fits"}, # 4 exts
    {'id': "4460981", 'filename':"1311160005010021.fits"} # 5 exts
])
def test_galah_guess(remote_data_path):
    spectra = SpectrumList.read(remote_data_path)
    # Should be main spectra, without sky, and normalised (not in 4 ext)
    nspec = len(spectra)
    if spectra[0].meta["galah_hdu_format"] == 4:
        assert nspec == 2
    elif spectra[0].meta["galah_hdu_format"] == 5:
        assert nspec == 3
    else:
        assert False, "Unknown format"

    # main spectra
    assert spectra[0].flux.unit == u.count
    assert spectra[0].spectral_axis.unit == u.Angstrom
    assert isinstance(spectra[0].uncertainty, StdDevUncertainty)
    assert spectra[0].meta.get("label") is not None
    assert spectra[0].meta.get("header") is not None

    # normalised
    if nspec == 3:
        assert spectra[1].flux.unit == u.Unit('') # dimensionless
        assert spectra[1].spectral_axis.unit == u.Angstrom
        assert spectra[1].uncertainty is None
        assert spectra[1].meta.get("label") == "normalised spectra"
        assert spectra[1].meta.get("header") is not None

    # No sky
    if nspec == 3:
        assert spectra[2].spectral_axis.unit == u.Angstrom
        assert spectra[2].flux.unit == u.count
        assert isinstance(spectra[2].uncertainty, StdDevUncertainty)
        assert spectra[2].meta.get("label") is not None
        assert spectra[2].meta.get("header") is not None
    else:
        assert spectra[1].spectral_axis.unit == u.Angstrom
        assert spectra[1].flux.unit == u.count
        assert isinstance(spectra[1].uncertainty, StdDevUncertainty)
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None


# We cannot use remote_access directly in the MIRI MRS tests, because
# the test functions are called once for every file in the remote access
# list, but we need to feed the entire set to the functions at once. We
# store the file names in a list and use a dummy test method to load the
# list via the remote_access machinery.
#
# MIRI MRS 1D data sets are comprised of 12 files.
# We add one bad file to test the skip/warn on missing file functionality
filename_list = ["bad_file.fits"]
@remote_access([
    {'id': '5082863', 'filename': 'combine_dithers_all_exposures_ch1-long_x1d.fits'},
    {'id': '5082863', 'filename': 'combine_dithers_all_exposures_ch1-medium_x1d.fits'},
    {'id': '5082863', 'filename': 'combine_dithers_all_exposures_ch1-short_x1d.fits'},
    {'id': '5082863', 'filename': 'combine_dithers_all_exposures_ch2-long_x1d.fits'},
    {'id': '5082863', 'filename': 'combine_dithers_all_exposures_ch2-medium_x1d.fits'},
    {'id': '5082863', 'filename': 'combine_dithers_all_exposures_ch2-short_x1d.fits'},
    {'id': '5082863', 'filename': 'combine_dithers_all_exposures_ch3-long_x1d.fits'},
    {'id': '5082863', 'filename': 'combine_dithers_all_exposures_ch3-medium_x1d.fits'},
    {'id': '5082863', 'filename': 'combine_dithers_all_exposures_ch3-short_x1d.fits'},
    {'id': '5082863', 'filename': 'combine_dithers_all_exposures_ch4-long_x1d.fits'},
    {'id': '5082863', 'filename': 'combine_dithers_all_exposures_ch4-medium_x1d.fits'},
    {'id': '5082863', 'filename': 'combine_dithers_all_exposures_ch4-short_x1d.fits'},
])
def test_loaddata_miri_mrs(remote_data_path):
    filename_list.append(remote_data_path)


# loading from a list of file names
@pytest.mark.remote_data
def test_spectrum_list_names_miri_mrs(caplog):

    # Format is explicitly set
    with pytest.raises(FileNotFoundError):
        specs = SpectrumList.read(filename_list, format="JWST x1d MIRI MRS")

    # Skip missing file silently
    specs = SpectrumList.read(filename_list, format="JWST x1d MIRI MRS", missing='silent')

    assert len(specs) == 12
    for spec in specs:
        assert isinstance(spec, Spectrum1D)
        assert spec.spectral_axis.unit == u.micron

    # Warn about missing file
    specs = SpectrumList.read(filename_list, format="JWST x1d MIRI MRS", missing='warn')

    assert "Failed to load bad_file.fits: FileNotFoundError" in caplog.text

    assert len(specs) == 12
    for spec in specs:
        assert isinstance(spec, Spectrum1D)
        assert spec.spectral_axis.unit == u.Unit("um")

    # Auto-detect format
    specs = SpectrumList.read(filename_list[1:])

    assert len(specs) == 12
    for spec in specs:
        assert isinstance(spec, Spectrum1D)
        assert spec.spectral_axis.unit == u.micron


# loading from a directory via glob
@pytest.mark.remote_data
def test_spectrum_list_directory_miri_mrs(tmpdir):

    # copy files to temp dir. We cannot use the directory generated by
    # remote_access, because it may have variable structure from run to
    # run. And also because temp directories created in previous runs
    # may still be hanging around. This precludes the use of commonpath()
    tmp_dir = tmpdir.strpath
    for file_path in filename_list[1:]:
        shutil.copy(file_path, tmp_dir)

    specs = SpectrumList.read(tmp_dir)

    assert len(specs) == 12
    for spec in specs:
        assert isinstance(spec, Spectrum1D)
        assert spec.spectral_axis.unit == u.micron


# x1d and c1d assorted files

@remote_access([
    {'id': "5394931", 'filename':"jw00623-c1012_t002_miri_p750l_x1d.fits"},                  # pipeline 1.2.3
    {'id': "5394931", 'filename':"jw00787-o014_s00002_niriss_f150w-gr150c-gr150r_c1d.fits"}, # pipeline 1.2.3
    {'id': "5394931", 'filename':"jw00623-o057_t008_miri_ch1-long_x1d.fits"},                # pipeline 1.3.1
    {'id': "5394931", 'filename':"jw00626-o064_t007_nirspec_g235h-f170lp_x1d.fits"},         # pipeline 1.3.1
])
def test_jwst_x1d_c1d(remote_data_path):

    data = Spectrum1D.read(remote_data_path)

    assert isinstance(data, Spectrum1D)
    assert data.shape in [(388,), (5,), (1091,), (3843,)]
    assert data.unit == u.Jy
    assert data.spectral_axis.unit == u.um


# utility functions to be used with list comprehension in SpectrumList checking
def assert_multi_isinstance(a, b):
    assert isinstance(a, b)

def assert_multi_equals(a, b):
    assert a == b


@remote_access([
    {'id': "5394931", 'filename':"jw00624-o027_s00001_nircam_f356w-grismr_x1d.fits"}, # pipeline 1.2.3
])
def test_jwst_nircam_x1d_multi_v1_2_3(remote_data_path):

    data = SpectrumList.read(remote_data_path)

    assert isinstance(data, SpectrumList)
    assert len(data) == 3
    [assert_multi_isinstance(d, Spectrum1D) for d in data]
    [assert_multi_equals(d.shape, r) for d,r in zip(data, [(459,), (336,), (962,)])]
    [assert_multi_equals(d.unit, u.Jy) for d in data]
    [assert_multi_equals(d.spectral_axis.unit, u.um) for d in data]


@remote_access([
    {'id': "5394931", 'filename':"jw00660-o016_s00002_nircam_f444w-grismr_x1d.fits"}, # pipeline 1.3.1
])
def test_jwst_nircam_x1d_multi_v1_3_1(remote_data_path):

    data = SpectrumList.read(remote_data_path)

    assert isinstance(data, SpectrumList)
    assert len(data) == 4
    [assert_multi_isinstance(d, Spectrum1D) for d in data]
    [assert_multi_equals(d.shape, r) for d,r in zip(data, [(1166,), (786,), (1157,), (795,)])]
    [assert_multi_equals(d.unit, u.Jy) for d in data]
    [assert_multi_equals(d.spectral_axis.unit, u.um) for d in data]


@remote_access([
    {'id': "5394931", 'filename':"jw00776-o003_s00083_nircam_f322w2-grismr_c1d.fits"}, # pipeline 1.2.3
])
def test_jwst_nircam_c1d_v1_2_3(remote_data_path):

    data = SpectrumList.read(remote_data_path)

    assert isinstance(data, SpectrumList)
    assert len(data) == 2
    [assert_multi_isinstance(d, Spectrum1D) for d in data]
    [assert_multi_equals(d.shape, r) for d,r in zip(data, [(133,), (1139,)])]
    [assert_multi_equals(d.unit, u.MJy/u.sr) for d in data]
    [assert_multi_equals(d.spectral_axis.unit, u.um) for d in data]


@remote_access([
    {'id': "5394931", 'filename':"jw00625-o018_s00001_niriss_f090w-gr150c_c1d.fits"}, # pipeline 1.2.3
])
def test_jwst_niriss_c1d_v1_2_3(remote_data_path):

    data = SpectrumList.read(remote_data_path)

    assert isinstance(data, SpectrumList)
    assert len(data) == 2
    [assert_multi_isinstance(d, Spectrum1D) for d in data]
    [assert_multi_equals(d.shape, r) for d,r in zip(data, [(56,), (107,)])]
    [assert_multi_equals(d.unit, u.Jy) for d in data]
    [assert_multi_equals(d.spectral_axis.unit, u.um) for d in data]
