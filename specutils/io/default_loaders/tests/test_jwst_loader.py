import numpy as np
from astropy.io import fits
from astropy.table import Table
import pytest

from specutils import Spectrum1D, SpectrumList


def create_spectrum_hdu(data_len, srctype=None):
    # Create a minimal header for the purposes of testing

    data = np.random.random((data_len, 5))
    table = Table(data=data, names=['WAVELENGTH', 'FLUX', 'ERROR', 'SURF_BRIGHT',
        'SB_ERROR'])

    hdu = fits.BinTableHDU(table, name='EXTRACT1D')
    hdu.header['TUNIT1'] = 'um'
    hdu.header['TUNIT2'] = 'Jy'
    hdu.header['TUNIT3'] = 'Jy'
    hdu.header['TUNIT4'] = 'MJy/sr'
    hdu.header['TUNIT5'] = 'MJy/sr'
    hdu.header['SRCTYPE'] = srctype

    return hdu


@pytest.fixture
def x1d_single():
    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())
    hdulist["PRIMARY"].header["TELESCOP"] = "JWST"
    # Add several BinTableHDUs that contain spectral data
    hdulist.append(create_spectrum_hdu(100, 'POINT'))
    # Mock the ASDF extension
    hdulist.append(fits.BinTableHDU(name='ASDF'))

    return hdulist


@pytest.fixture
def x1d_multi():
    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())
    hdulist["PRIMARY"].header["TELESCOP"] = "JWST"
    # Add several BinTableHDUs that contain spectral data
    hdulist.append(create_spectrum_hdu(100, 'POINT'))
    hdulist.append(create_spectrum_hdu(120, 'EXTENDED'))
    hdulist.append(create_spectrum_hdu(110, 'POINT'))
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


def test_jwst_x1d_singel_loader_fail_on_multi(tmpdir, x1d_multi):
    """Make sure Spectrum1D.read on JWST x1d with many spectra errors out"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    x1d_multi.writeto(tmpfile)

    with pytest.raises(RuntimeError, match="SpectrumList"):
        Spectrum1D.read(tmpfile, format='JWST x1d')


@pytest.mark.parametrize("srctype", [None, "UNKNOWN"])
def test_jwst_loader_fail(tmpdir, x1d_single, srctype):
    """Check that the loader fails when SRCTYPE is not set or is UNKNOWN"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    hdulist = x1d_single
    # Add a spectrum with unknown SRCTYPE
    hdulist.append(create_spectrum_hdu(100, srctype))
    hdulist.writeto(tmpfile)

    with pytest.raises(RuntimeError, match="^Keyword"):
        SpectrumList.read(tmpfile, format='JWST x1d')
