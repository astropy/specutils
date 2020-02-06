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


def test_jwst_loader(tmpdir):

    tmpfile = str(tmpdir.join('jwst.fits'))

    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())
    hdulist["PRIMARY"].header["TELESCOP"] = "JWST"
    # Add several BinTableHDUs that contain spectral data
    hdulist.append(create_spectrum_hdu(100, 'POINT'))
    hdulist.append(create_spectrum_hdu(120, 'EXTENDED'))
    hdulist.append(create_spectrum_hdu(110, 'POINT'))
    # Mock the ASDF extension
    hdulist.append(fits.BinTableHDU(name='ASDF'))
    hdulist.writeto(tmpfile)

    data = SpectrumList.read(tmpfile, format='JWST')
    assert len(data) == 3

    for item in data:
        assert isinstance(item, Spectrum1D)

    assert data[0].shape == (100,)
    assert data[1].shape == (120,)
    assert data[2].shape == (110,)

def test_jwst_loader_fail(tmpdir):

    tmpfile = str(tmpdir.join('jwst.fits'))

    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())
    hdulist["PRIMARY"].header["TELESCOP"] = "JWST"
    # Add several BinTableHDUs that contain spectral data
    hdulist.append(create_spectrum_hdu(100, 'UNKNOWN'))
    # Mock the ASDF extension
    hdulist.append(fits.BinTableHDU(name='ASDF'))
    hdulist.writeto(tmpfile)

    with pytest.raises(RuntimeError, match="^Keyword"):
        SpectrumList.read(tmpfile, format='JWST')
