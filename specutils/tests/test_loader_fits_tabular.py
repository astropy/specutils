import os
import warnings

import pytest
import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.wcs import FITSFixedWarning
from astropy.nddata import StdDevUncertainty
from astropy.tests.helper import quantity_allclose

from numpy.testing import assert_allclose

from .. import Spectrum1D


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


@pytest.mark.parametrize("spectral_axis",
                         ['wavelength', 'frequency', 'energy', 'wavenumber'])
def test_tabular_fits_writer(tmp_path, spectral_axis):
    wlu = {'wavelength': u.AA, 'frequency': u.GHz, 'energy': u.eV,
           'wavenumber': u.cm**-1}
    # Create a small data set
    disp = np.arange(1, 1.1, 0.01) * wlu[spectral_axis]
    flux = np.ones(len(disp)) * 1.e-14 * u.Jy
    unc = StdDevUncertainty(0.01 * flux)
    if spectral_axis not in ('wavelength', ):
        disp = np.flip(disp)

    spectrum = Spectrum1D(flux=flux, spectral_axis=disp, uncertainty=unc)
    tmpfile = str(tmp_path / '_tst.fits')
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
def test_tabular_fits_multid(tmp_path, ndim, spectral_axis):
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
    tmpfile = str(tmp_path / '_tst.fits')
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


def test_tabular_fits_header(tmp_path):
    # Create a small data set + header with reserved FITS keywords
    disp = np.linspace(1, 1.2, 21) * u.AA
    flux = np.random.normal(0., 1.0e-14, disp.shape[0]) * u.Jy
    hdr = fits.header.Header({'TELESCOP': 'Leviathan', 'APERTURE': 1.8,
                              'OBSERVER': 'Parsons'})

    spectrum = Spectrum1D(flux=flux, spectral_axis=disp, meta={'header': hdr})
    tmpfile = str(tmp_path / '_tst.fits')
    spectrum.write(tmpfile, format='tabular-fits')

    print(f"spetrum meta: \n {repr(spectrum.meta['header'])}")
    # Read it in and check against the original
    with fits.open(tmpfile) as hdulist:
        print(f'\nhdu info:\n')
        hdulist.info()
        print(f'\nhdu0 header:\n{repr(hdulist[0].header)}')
        print(f'\nhdu1 header:\n{repr(hdulist[1].header)}')

        assert hdulist[0].header['OBSERVER'] == 'Parsons'
    
        # keys relevant to datashape are in HDU 1
        assert hdulist[1].header['NAXIS'] == 2
        assert hdulist[1].header['NAXIS2'] == disp.shape[0]


def test_tabular_fits_update_header(tmp_path):
    disp = np.linspace(1, 1.2, 21) * u.AA
    flux = np.random.normal(0., 1.0e-14, disp.shape[0]) * u.erg / (u.s * u.cm**2 * u.AA)
    hdr = fits.header.Header({'TELESCOP': 'Crystal', 'OBSERVER': 'Cruz'})
    spec = Spectrum1D(flux=flux, spectral_axis=disp, meta={'header': hdr})
    tmpfile = str(tmp_path / '_tst2.fits')
    spec.write(tmpfile, format='tabular-fits')
    spectrum = Spectrum1D.read(tmpfile, format='tabular-fits')

    # Now write with updated header information from spectrum.meta
    spectrum.meta.update({'OBSERVER': 'Rosse', 'EXPTIME': 32.1, 'NAXIS2': 12})
    spectrum.write(tmpfile, format='tabular-fits', overwrite=True,
                   update_header=True)

    with fits.open(tmpfile) as hdulist:
        assert hdulist[0].header['OBSERVER'] == 'Rosse'
        #  assert hdulist[1].header['NAXIS2'] == disp.shape[0]
        assert_allclose(hdulist[1].header['EXPTIME'], 3.21e1)

    # Test that unsupported types (dict) are not added to written header
    spectrum.meta['MYHEADER'] = {'OBSDATE': '1848-02-26', 'TARGET': 'M51'}
    spectrum.write(tmpfile, format='tabular-fits', overwrite=True,
                   update_header=True)

    with fits.open(tmpfile) as hdulist:
        assert 'MYHEADER' not in hdulist[0].header
        assert 'MYHEADER' not in hdulist[1].header
        assert 'OBSDATE' not in hdulist[0].header
        assert 'OBSDATE' not in hdulist[1].header


@pytest.mark.filterwarnings("ignore:The unit 'Angstrom' has been deprecated")
def test_tabular_fits_autowrite(tmp_path):
    """Test writing of Spectrum1D with automatic selection of BINTABLE format."""
    disp = np.linspace(1, 1.2, 21) * u.AA
    flux = np.random.normal(0., 1.0e-14, disp.shape[0]) * u.W / (u.m**2 * u.AA)
    hdr = fits.header.Header({'TELESCOP': 'Leviathan', 'APERTURE': 1.8,
                              'OBSERVER': 'Parsons'})

    spectrum = Spectrum1D(flux=flux, spectral_axis=disp, meta={'header': hdr})
    tmpfile = str(tmp_path / '_tst.fits')
    spectrum.write(tmpfile)

    # Read it in and check against the original
    with fits.open(tmpfile) as hdulist:
        assert hdulist[1].header['NAXIS'] == 2
        assert hdulist[1].header['NAXIS2'] == disp.shape[0]

    # Trigger exception for illegal HDU (primary HDU only accepts IMAGE_HDU)
    with pytest.raises(ValueError, match=r'FITS does not support BINTABLE'):
        spectrum.write(tmpfile, format='tabular-fits', overwrite=True, hdu=0)

    # Test automatic selection of wcs1d format, which will fail without suitable wcs
    with pytest.raises(ValueError, match=r'Only Spectrum1D objects with valid WCS'):
        spectrum.write(tmpfile, overwrite=True, hdu=0)

    tmpfile = str(tmp_path / '_wcs.fits')
    with pytest.raises(ValueError, match=r'Only Spectrum1D objects with valid WCS'):
        spectrum.write(tmpfile, overwrite=True)


@pytest.mark.skipif('sys.platform.startswith("win")',
                    reason='Uncertain availability of compression utilities')
@pytest.mark.parametrize('compress', ['gzip', 'bzip2', 'xz'])
def test_tabular_fits_compressed(compress, tmp_path):
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
    tmpfile = str(tmp_path / '_tst.fits')
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
