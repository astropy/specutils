import numpy as np
from astropy.io import fits
from astropy.table import Table

from specutils import Spectrum1D, SpectrumList


def create_spectrum_hdu(data_len):
    # Create a minimal header for the purposes of testing

    data = np.random.random((data_len, 3))
    table = Table(data=data, names=['WAVELENGTH', 'FLUX', 'ERROR'])

    hdu = fits.BinTableHDU(table, name='EXTRACT1D')
    hdu.header['TUNIT1'] = 'um'
    hdu.header['TUNIT2'] = 'mJy'
    hdu.header['TUNIT3'] = 'mJy'

    return hdu


def test_jwst_loader(tmpdir):

    tmpfile = str(tmpdir.join('jwst.fits'))

    hdulist = fits.HDUList()
    # Make sure the file has a primary HDU
    hdulist.append(fits.PrimaryHDU())
    # Add several BinTableHDUs that contain spectral data
    hdulist.append(create_spectrum_hdu(100))
    hdulist.append(create_spectrum_hdu(120))
    hdulist.append(create_spectrum_hdu(110))
    # JWST data product will always contain an ASDF header which is a BinTable
    hdulist.append(fits.BinTableHDU(name='ASDF'))
    hdulist.writeto(tmpfile)

    data = SpectrumList.read(tmpfile, format='JWST')
    assert len(data) == 3

    for item in data:
        assert isinstance(item, Spectrum1D)

    assert data[0].shape == (100,)
    assert data[1].shape == (120,)
    assert data[2].shape == (110,)
