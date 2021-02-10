import numpy as np
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from astropy.io.registry import IORegistryError
from astropy.utils.exceptions import AstropyUserWarning
from astropy.modeling import models
from astropy import coordinates as coord
import gwcs.coordinate_frames as cf
from gwcs.wcs import WCS
import pytest

from specutils import Spectrum1D, SpectrumList


# The x1d reader tests --------------------------

def create_spectrum_hdu(data_len, srctype=None, ver=1):
    """Mock a JWST x1d BinTableHDU"""
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
    hdu.ver = ver

    return hdu


@pytest.fixture(scope="function")
def x1d_single():
    """Mock a JWST x1d HDUList with a single spectrum"""
    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())
    hdulist["PRIMARY"].header["TELESCOP"] = ("JWST", "comment")
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


def test_jwst_x1d_multi_reader(tmpdir, x1d_multi):
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


def test_jwst_x1d_single_reader(tmpdir, x1d_single):
    """Test Spectrum1D.read for JWST x1d data"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    x1d_single.writeto(tmpfile)

    data = Spectrum1D.read(tmpfile, format='JWST x1d')
    assert type(data) is Spectrum1D
    assert data.shape == (100,)


def test_jwst_x1d_single_reader_no_format(tmpdir, x1d_single):
    """Test Spectrum1D.read for JWST x1d data without format arg"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    x1d_single.writeto(tmpfile)

    data = Spectrum1D.read(tmpfile)
    assert type(data) is Spectrum1D
    assert data.shape == (100,)
    assert data.unit == u.Jy
    assert data.spectral_axis.unit == u.um


def test_jwst_x1d_multi_reader_no_format(tmpdir, x1d_multi):
    """Test Spectrum1D.read for JWST x1d data without format arg"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    x1d_multi.writeto(tmpfile)

    data = SpectrumList.read(tmpfile)
    assert type(data) is SpectrumList
    assert len(data) == 3

    for item in data:
        assert isinstance(item, Spectrum1D)


def test_jwst_x1d_multi_reader_check_units(tmpdir, x1d_multi):
    """Test units for Spectrum1D.read for JWST x1d data"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    x1d_multi.writeto(tmpfile)

    data = SpectrumList.read(tmpfile)
    assert data[0].unit == u.Jy
    assert data[1].unit == u.MJy / u.sr
    assert data[2].unit == u.Jy


def test_jwst_x1d_reader_meta(tmpdir, x1d_single):
    """Test that the Primary and EXTRACT1D extension headers are merged in meta"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    x1d_single.writeto(tmpfile)

    data = Spectrum1D.read(tmpfile)
    assert ('TELESCOP', 'JWST') in data.meta['header'].items()
    assert ('SRCTYPE', 'POINT') in data.meta['header'].items()


def test_jwst_x1d_single_reader_fail_on_multi(tmpdir, x1d_multi):
    """Make sure Spectrum1D.read on JWST x1d with many spectra errors out"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    x1d_multi.writeto(tmpfile)

    with pytest.raises(IORegistryError):
        Spectrum1D.read(tmpfile)


@pytest.mark.parametrize("srctype", [None, "UNKNOWN"])
def test_jwst_reader_fail(tmpdir, x1d_single, srctype):
    """Check that the reader fails when SRCTYPE is not set or is UNKNOWN"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    hdulist = x1d_single
    # Add a spectrum with bad SRCTYPE (mutate the fixture)
    hdulist.append(create_spectrum_hdu(100, srctype, ver=2))
    hdulist.writeto(tmpfile)

    with pytest.raises(RuntimeError, match="^Keyword"):
        SpectrumList.read(tmpfile, format='JWST x1d multi')


@pytest.mark.xfail(reason="JWST loader no longer attempts to auto-find flux column.")
def test_jwst_reader_warning_stddev(tmpdir, x1d_single):
    """Check that the reader raises warning when stddev is zeros"""
    tmpfile = str(tmpdir.join('jwst.fits'))
    hdulist = x1d_single
    # Put zeros in ERROR column
    hdulist["EXTRACT1D"].data["ERROR"] = 0
    hdulist.writeto(tmpfile)

    with pytest.warns(Warning) as record:
        Spectrum1D.read(tmpfile)
        for r in record:
            if r.message is AstropyUserWarning:
                assert "Standard Deviation has values of 0" in r.message


# The s2d/s3d reader tests -------------------------------

@pytest.fixture
def generate_wcs_transform():
    def _generate_wcs_transform(dispaxis):
        """Create mock gwcs.WCS object for resampled s2d data"""
        detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
        icrs = cf.CelestialFrame(name='icrs', reference_frame=coord.ICRS(),
            axes_order=(0, 1), unit=(u.deg, u.deg), axes_names=('RA', 'DEC'))
        spec = cf.SpectralFrame(name='spec', axes_order=(2,), unit=(u.micron,),
            axes_names=('lambda',))
        world = cf.CompositeFrame(name="world", frames=[icrs, spec])

        if dispaxis == 1:
            mapping = models.Mapping((0, 1, 0))
        if dispaxis == 2:
            mapping = models.Mapping((0, 1, 1))

        transform = mapping | (models.Const1D(42) & models.Const1D(42)
            & (models.Shift(30) | models.Scale(0.1)))
        pipeline = [(detector, transform),
                    (world, None)]
        wcs = WCS(pipeline)

        return wcs

    return _generate_wcs_transform


@pytest.fixture
def s2d_single(generate_wcs_transform):
    pytest.importorskip("jwst")
    from jwst.datamodels import MultiSlitModel, SlitModel
    from jwst.assign_wcs.util import wcs_bbox_from_shape

    shape = (10, 100)
    dispaxis = 1

    model = MultiSlitModel()
    sm = SlitModel(shape)
    sm.data
    model.slits.append(sm)
    for slit in model.slits:
        slit.meta.wcs = generate_wcs_transform(dispaxis)
        slit.meta.wcs.bounding_box = wcs_bbox_from_shape(shape)
        slit.meta.wcsinfo.dispersion_direction = dispaxis

    model.meta.telescope = "JWST"

    return model


@pytest.fixture(params=[(5, 100), (100, 8)])
def s2d_multi(generate_wcs_transform, request):
    pytest.importorskip("jwst")
    from jwst.datamodels import SlitModel, MultiSlitModel
    from jwst.assign_wcs.util import wcs_bbox_from_shape

    shape = request.param
    if shape[0] < shape[1]:
        dispaxis = 1
    else:
        dispaxis = 2

    model = MultiSlitModel()
    sm = SlitModel(shape)
    sm.data
    model.slits.append(sm)
    model.slits.append(sm)
    for slit in model.slits:
        slit.meta.wcs = generate_wcs_transform(dispaxis)
        slit.meta.wcs.bounding_box = wcs_bbox_from_shape(shape)
        slit.meta.wcsinfo.dispersion_direction = dispaxis
        slit.meta.bunit_data = "Jy"
        slit.meta.bunit_err = "Jy"

    return model


@pytest.mark.xfail(reason="Needs investigation! See #717")
def test_jwst_s2d_reader(tmpdir, s2d_single):
    path = str(tmpdir.join("test.fits"))
    model = s2d_single
    model.save(path)

    spec = Spectrum1D.read(path)
    assert hasattr(spec, "spectral_axis")
    assert spec.unit == u.dimensionless_unscaled


def test_jwst_s2d_multi_reader(tmpdir, s2d_multi):
    path = str(tmpdir.join("test.fits"))
    model = s2d_multi
    model.save(path)

    speclist = SpectrumList.read(path, format="JWST s2d multi")
    assert len(speclist) == 2
    assert hasattr(speclist[0], "spectral_axis")
    assert speclist[1].unit == u.Jy
