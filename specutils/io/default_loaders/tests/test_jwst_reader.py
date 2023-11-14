import gwcs.coordinate_frames as cf
import numpy as np
import pytest
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import fits
from astropy.io.registry import IORegistryError
from astropy.modeling import models
from astropy.table import Table
from astropy.utils.exceptions import AstropyUserWarning
from gwcs.wcs import WCS

from specutils import Spectrum1D, SpectrumList

try:
    from stdatamodels import asdf_in_fits
except ImportError:
    HAS_STDATAMODELS = False
else:
    HAS_STDATAMODELS = True


# The c1d/x1d reader tests --------------------------

def create_spectrum_hdu(data_len, srctype=None, ver=1, name='EXTRACT1D'):
    """Mock a JWST x1d BinTableHDU"""
    np.random.seed(20)
    data = np.random.random((data_len, 5))

    # make sure spectral axis is sorted
    data = data[data[:, 0].argsort()]

    table = Table(data=data, names=['WAVELENGTH', 'FLUX', 'ERROR', 'SURF_BRIGHT',
        'SB_ERROR'])

    hdu = fits.BinTableHDU(table, name=name)
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
def spec_single(request):
    """Mock a JWST c1d/x1d HDUList with a single spectrum"""
    name = request.param
    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())
    hdulist["PRIMARY"].header["TELESCOP"] = ("JWST", "comment")
    # Add a BinTableHDU that contains spectral data
    hdulist.append(create_spectrum_hdu(100, 'POINT', ver=1, name=name))
    # Mock the ASDF extension
    hdulist.append(fits.BinTableHDU(name='ASDF'))

    return hdulist


@pytest.fixture(scope="function")
def spec_multi(request):
    """Mock a JWST c1d/x1d multispec HDUList with 3 spectra"""
    name = request.param
    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())
    hdulist["PRIMARY"].header["TELESCOP"] = "JWST"
    # Add a few BinTableHDUs that contain spectral data
    hdulist.append(create_spectrum_hdu(100, 'POINT', ver=1, name=name))
    hdulist.append(create_spectrum_hdu(120, 'EXTENDED', ver=2, name=name))
    hdulist.append(create_spectrum_hdu(110, 'POINT', ver=3, name=name))
    # Mock the ASDF extension
    hdulist.append(fits.BinTableHDU(name='ASDF'))

    return hdulist


@pytest.mark.parametrize('spec_multi, format',
                         [('EXTRACT1D', 'JWST x1d multi'),
                          ('COMBINE1D', 'JWST c1d multi')], indirect=['spec_multi'])
def test_jwst_1d_multi_reader(tmp_path, spec_multi, format):
    """Test SpectrumList.read for JWST c1d/x1d multi data"""
    tmpfile = str(tmp_path / 'jwst.fits')
    spec_multi.writeto(tmpfile)

    data = SpectrumList.read(tmpfile, format=format)
    assert type(data) is SpectrumList
    assert len(data) == 3

    for item in data:
        assert isinstance(item, Spectrum1D)

    assert data[0].shape == (100,)
    assert data[1].shape == (120,)
    assert data[2].shape == (110,)


@pytest.mark.parametrize('spec_single, format',
                         [('EXTRACT1D', 'JWST x1d'),
                          ('COMBINE1D', 'JWST c1d')], indirect=['spec_single'])
def test_jwst_1d_single_reader(tmp_path, spec_single, format):
    """Test Spectrum1D.read for JWST x1d data"""
    tmpfile = str(tmp_path / 'jwst.fits')
    spec_single.writeto(tmpfile)

    data = Spectrum1D.read(tmpfile, format=format)
    assert type(data) is Spectrum1D
    assert data.shape == (100,)


@pytest.mark.parametrize("srctype", [None, "UNKNOWN"])
def test_jwst_srctpye_defaults(tmp_path, x1d_single, srctype):
    """ Test """
    tmpfile = str(tmp_path / 'jwst.fits')

    # Add a spectrum with missing or UNKNOWN SRCTYPE (mutate the fixture)
    x1d_single['EXTRACT1D'].header['SRCTYPE'] == srctype
    x1d_single.writeto(tmpfile)

    data = Spectrum1D.read(tmpfile, format='JWST x1d')
    assert type(data) is Spectrum1D
    assert data.shape == (100,)
    assert x1d_single['EXTRACT1D'].header['SRCTYPE'] == "POINT"


@pytest.mark.parametrize('spec_single', ['EXTRACT1D', 'COMBINE1D'], indirect=['spec_single'])
def test_jwst_1d_single_reader_no_format(tmp_path, spec_single):
    """Test Spectrum1D.read for JWST c1d/x1d data without format arg"""
    tmpfile = str(tmp_path / 'jwst.fits')
    spec_single.writeto(tmpfile)

    data = Spectrum1D.read(tmpfile)
    assert type(data) is Spectrum1D
    assert data.shape == (100,)
    assert data.unit == u.Jy
    assert data.spectral_axis.unit == u.um


@pytest.mark.parametrize('spec_multi', ['EXTRACT1D', 'COMBINE1D'], indirect=['spec_multi'])
def test_jwst_1d_multi_reader_no_format(tmp_path, spec_multi):
    """Test Spectrum1D.read for JWST c1d/x1d data without format arg"""
    tmpfile = str(tmp_path / 'jwst.fits')
    spec_multi.writeto(tmpfile)

    data = SpectrumList.read(tmpfile)
    assert type(data) is SpectrumList
    assert len(data) == 3

    for item in data:
        assert isinstance(item, Spectrum1D)


@pytest.mark.parametrize('spec_multi', ['EXTRACT1D', 'COMBINE1D'], indirect=['spec_multi'])
def test_jwst_1d_multi_reader_check_units(tmp_path, spec_multi):
    """Test units for Spectrum1D.read for JWST c1d/x1d data"""
    tmpfile = str(tmp_path / 'jwst.fits')
    spec_multi.writeto(tmpfile)

    data = SpectrumList.read(tmpfile)
    assert data[0].unit == u.Jy
    assert data[1].unit == u.MJy / u.sr
    assert data[2].unit == u.Jy


@pytest.mark.parametrize('spec_single', ['EXTRACT1D', 'COMBINE1D'], indirect=['spec_single'])
def test_jwst_1d_reader_meta(tmp_path, spec_single):
    """Test that the Primary and COMBINE1D/EXTRACT1D extension headers are merged in meta"""
    tmpfile = str(tmp_path / 'jwst.fits')
    spec_single.writeto(tmpfile)

    data = Spectrum1D.read(tmpfile)
    assert ('TELESCOP', 'JWST') in data.meta['header'].items()
    assert ('SRCTYPE', 'POINT') in data.meta['header'].items()


@pytest.mark.parametrize('spec_multi', ['EXTRACT1D', 'COMBINE1D'], indirect=['spec_multi'])
def test_jwst_1d_single_reader_fail_on_multi(tmp_path, spec_multi):
    """Make sure Spectrum1D.read on JWST c1d/x1d with many spectra errors out"""
    tmpfile = str(tmp_path / 'jwst.fits')
    spec_multi.writeto(tmpfile)

    with pytest.raises(IORegistryError):
        Spectrum1D.read(tmpfile)


@pytest.mark.parametrize("srctype", ["BADVAL"])
def test_jwst_reader_fail(tmp_path, x1d_single, srctype):
    """Check that the reader fails when SRCTYPE is a BADVAL"""
    tmpfile = str(tmp_path / 'jwst.fits')
    hdulist = x1d_single
    # Add a spectrum with bad SRCTYPE (mutate the fixture)
    hdulist.append(create_spectrum_hdu(100, srctype, ver=2))
    hdulist.writeto(tmpfile)

    with pytest.raises(RuntimeError, match="^Keyword"):
        SpectrumList.read(tmpfile, format='JWST x1d multi')


@pytest.mark.xfail(reason="JWST loader no longer attempts to auto-find flux column.")
def test_jwst_reader_warning_stddev(tmp_path, x1d_single):
    """Check that the reader raises warning when stddev is zeros"""
    tmpfile = str(tmp_path / 'jwst.fits')
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
def test_jwst_s2d_reader(tmp_path, s2d_single):
    path = str(tmp_path / "test.fits")
    model = s2d_single
    model.save(path)

    spec = Spectrum1D.read(path)
    assert hasattr(spec, "spectral_axis")
    assert spec.unit == u.dimensionless_unscaled


def test_jwst_s2d_multi_reader(tmp_path, s2d_multi):
    path = str(tmp_path / "test.fits")
    model = s2d_multi
    model.save(path)

    speclist = SpectrumList.read(path, format="JWST s2d multi")
    assert len(speclist) == 2
    assert hasattr(speclist[0], "spectral_axis")
    assert speclist[1].unit == u.Jy


# The s3d reader tests -------------------------------

def generate_s3d_wcs():
    """ create a fake gwcs for a cube """
    # create input /output frames
    detector = cf.CoordinateFrame(name='detector', axes_order=(0,1,2), axes_names=['x', 'y', 'z'],
                                  axes_type=['spatial', 'spatial', 'spatial'], naxes=3,
                                  unit=['pix', 'pix', 'pix'])
    sky = cf.CelestialFrame(reference_frame=coord.ICRS(), name='sky', axes_names=("RA", "DEC"))
    spec = cf.SpectralFrame(name='spectral', unit=['um'], axes_names=['wavelength'], axes_order=(2,))
    world = cf.CompositeFrame(name="world", frames=[sky, spec])

    # create fake transform to at least get a bounding box
    # for the s3d jwst loader

    # shape 30,10,10 (spec, y, x)
    crpix1, crpix2, crpix3 = 5, 5, 15  # (x, y, spec)
    crval1, crval2, crval3 = 1, 1, 1
    cdelt1, cdelt2, cdelt3 = 0.01, 0.01, 0.05

    shift = models.Shift(-crpix2) & models.Shift(-crpix1)
    scale = models.Multiply(cdelt2) & models.Multiply(cdelt1)
    proj = models.Pix2Sky_TAN()
    skyrot = models.RotateNative2Celestial(crval2, 90 + crval1, 180)
    celestial = shift | scale | proj | skyrot
    wave_model = models.Shift(-crpix3) | models.Multiply(cdelt3) | models.Shift(crval3)
    transform = models.Mapping((2, 0, 1)) | celestial & wave_model | models.Mapping((1, 2, 0))
    # bounding box based on shape (30,10,10) in test
    transform.bounding_box = ((0, 29), (0, 9), (0, 9))

    # create final wcs
    pipeline = [(detector, transform),
                (world, None)]
    return WCS(pipeline)


@pytest.fixture()
def tmp_asdf():
    # Create some data
    sequence = np.arange(100)
    squares  = sequence**2
    random = np.random.random(100)

    # Store the data in an arbitrarily nested dictionary
    tree = {
        'foo': 42,
        'name': 'Monty',
        'sequence': sequence,
        'powers': { 'squares' : squares },
        'random': random,
        'meta': {
            'wcs' : generate_s3d_wcs()
        }
    }

    yield tree
    tree = {}


def create_image_hdu(name='SCI', data=None, shape=None, hdrs=[], ndim=3):
    """ Mock an Image HDU """
    if data is None:
        if not shape:
            shape = [4, 2, 3] if ndim == 3 else [2, 2] if ndim == 2 else [2]
        data = np.zeros(shape)
    hdu = fits.ImageHDU(name=name, data=data, header=fits.Header(hdrs))
    hdu.ver = 1
    return hdu


@pytest.fixture(scope='function')
def cube(tmp_path, tmp_asdf):
    """ Mock a JWST s3d cube """
    prihdu = fits.PrimaryHDU()
    prihdu.header["TELESCOP"] = ("JWST", "comment")
    prihdu.header["FLUXEXT"] = ("ERR", "comment")
    prihdu.header["ERREXT"] = ("ERR", "comment")
    prihdu.header["MASKEXT"] = ("DQ", "comment")
    hdulist = fits.HDUList([prihdu])

    # Add ImageHDU for cubes
    shape = (30, 10, 10)
    hdulist.append(create_image_hdu(name='SCI', shape=shape, hdrs=[("BUNIT", 'MJy')]))
    hdulist.append(create_image_hdu(name='ERR', shape=shape,
                                    hdrs=[("BUNIT", 'MJy'), ('ERRTYPE', 'ERR')]))
    hdulist.append(create_image_hdu(name='DQ', shape=shape))

    # Mock the ASDF extension
    hdulist.append(fits.BinTableHDU(name='ASDF'))

    if HAS_STDATAMODELS:
        tmpfile = str(tmp_path / 'jwst_embedded_asdf.fits')
        asdf_in_fits.write(tmpfile, tmp_asdf, hdulist=hdulist, overwrite=True)

    return hdulist


@pytest.mark.skipif(not HAS_STDATAMODELS, reason="requires stdatamodels")
def test_jwst_s3d_single(tmp_path, cube):
    """Test Spectrum1D.read for JWST x1d data"""
    tmpfile = str(tmp_path / 'jwst_s3d.fits')
    cube.writeto(tmpfile)

    data = Spectrum1D.read(tmpfile, format='JWST s3d')
    assert type(data) is Spectrum1D
    assert data.shape == (10, 10, 30)
    assert data.uncertainty is not None
    assert data.mask is not None
    assert data.uncertainty.unit == 'MJy'
