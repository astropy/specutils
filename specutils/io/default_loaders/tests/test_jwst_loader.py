import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.io.registry import IORegistryError
from astropy.utils.exceptions import AstropyUserWarning
import pytest

from specutils import Spectrum1D, SpectrumList


def create_spectrum_hdu(data_len, srctype=None, ver=1):
    """Mock a JWST x1d BinTableHDU"""
    data = np.random.random((data_len, 5))
    table = Table(data=data,names=['WAVELENGTH', 'FLUX', 'ERROR', 'SURF_BRIGHT',
        'SB_ERROR'])

    hdu = fits.BinTableHDU(table, name='EXTRACT1D')
    hdu.header['TUNIT1'] = 'um'
    hdu.header['TUNIT2'] = 'Jy'
    hdu.header['TUNIT3'] = 'Jy'
    hdu.header['TUNIT4'] = 'MJy/sr'
    hdu.header['TUNIT5'] = 'MJy/sr'
    hdu.header['SRCTYPE'] = srctype
    hdu.ver = ver

    return hdu


@pytest.fixture(scope="function")
def x1d_single():
    """Mock a JWST x1d HDUList with a single spectrum"""
    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())
    hdulist["PRIMARY"].header["TELESCOP"] = "JWST"
    # Add a BinTableHDU that contains spectral data
    hdulist.append(create_spectrum_hdu(100, 'POINT', ver=1))
    # Mock the ASDF extension
    hdulist.append(fits.BinTableHDU(name='ASDF'))

    return hdulist


@pytest.fixture(scope="function")
def x1d_multi():
    """Mock a JWST x1d multispec HDUList with 3 spectra"""
    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())
    hdulist["PRIMARY"].header["TELESCOP"] = "JWST"
    # Add a few BinTableHDUs that contain spectral data
    hdulist.append(create_spectrum_hdu(100, 'POINT', ver=1))
    hdulist.append(create_spectrum_hdu(120, 'EXTENDED', ver=2))
    hdulist.append(create_spectrum_hdu(110, 'POINT', ver=3))
    # Mock the ASDF extension
    hdulist.append(fits.BinTableHDU(name='ASDF'))

    return hdulist


def test_jwst_x1d_multi_loader(tmpdir, x1d_multi):
    """Test SpectrumList.read for JWST x1d multi data"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    x1d_multi.writeto(tmpfile)

    data = SpectrumList.read(tmpfile, format='JWST x1d multi')
    assert type(data) is SpectrumList
    assert len(data) == 3

    for item in data:
        assert isinstance(item, Spectrum1D)

    assert data[0].shape == (100,)
    assert data[1].shape == (120,)
    assert data[2].shape == (110,)


def test_jwst_x1d_single_loader(tmpdir, x1d_single):
    """Test Spectrum1D.read for JWST x1d data"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    x1d_single.writeto(tmpfile)

    data = Spectrum1D.read(tmpfile, format='JWST x1d')
    assert type(data) is Spectrum1D
    assert data.shape == (100,)


def test_jwst_x1d_single_loader_no_format(tmpdir, x1d_single):
    """Test Spectrum1D.read for JWST x1d data without format arg"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    x1d_single.writeto(tmpfile)

    data = Spectrum1D.read(tmpfile)
    assert type(data) is Spectrum1D
    assert data.shape == (100,)


def test_jwst_x1d_multi_loader_no_format(tmpdir, x1d_multi):
    """Test Spectrum1D.read for JWST x1d data without format arg"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    x1d_multi.writeto(tmpfile)

    data = SpectrumList.read(tmpfile)
    assert type(data) is SpectrumList
    assert len(data) == 3

    for item in data:
        assert isinstance(item, Spectrum1D)

    assert data[0].shape == (100,)
    assert data[1].shape == (120,)
    assert data[2].shape == (110,)


def test_jwst_x1d_single_loader_fail_on_multi(tmpdir, x1d_multi):
    """Make sure Spectrum1D.read on JWST x1d with many spectra errors out"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    x1d_multi.writeto(tmpfile)

    with pytest.raises(IORegistryError):
        Spectrum1D.read(tmpfile)


@pytest.mark.parametrize("srctype", [None, "UNKNOWN"])
def test_jwst_loader_fail(tmpdir, x1d_single, srctype):
    """Check that the loader fails when SRCTYPE is not set or is UNKNOWN"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    hdulist = x1d_single
    # Add a spectrum with bad SRCTYPE (mutate the fixture)
    hdulist.append(create_spectrum_hdu(100, srctype, ver=2))
    hdulist.writeto(tmpfile)

    with pytest.raises(RuntimeError, match="^Keyword"):
        SpectrumList.read(tmpfile, format='JWST x1d multi')


def test_jwst_loader_warning_stddev(tmpdir, x1d_single):
    """Check that the loader raises warning when stddev is zeros"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    hdulist = x1d_single
    # Put zeros in ERROR column
    hdulist["EXTRACT1D"].data["ERROR"] = 0
    hdulist.writeto(tmpfile)

    with pytest.warns(AstropyUserWarning) as record:
        Spectrum1D.read(tmpfile)
        for r in record:
            if r.message is AstropyUserWarning:
                assert "Standard Deviation has values of 0" in r.message
