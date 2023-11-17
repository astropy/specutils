import gwcs.coordinate_frames as cf
import numpy as np
import pytest

from astropy.io import fits
from astropy.units import Unit, Quantity, Angstrom

from specutils import Spectrum1D, SpectrumList


def mwm_HDUList(n_spectra, format="visit"):
    """Mock a MWM HDUList"""
    np.random.seed(20)

    # TODO: i gave up on making this one.
    if format == "visit":
        hdulist = fits.open("./mwmVisit-0.5.0-70350000.fits")
    elif format == "star":
        hdulist = fits.open("./mwmStar-0.5.0-103020000.fits")

    return hdulist


def apStar_HDUList(n_spectra):
    """Mock an apStar HDUList of n_spectra spectra."""
    np.random.seed(20)

    # init primary hdu header
    hdr = fits.Header()
    hdr["FOOBAR"] = "barfoo"
    hdr["V_APRED"] = "q"
    hdr["APRED"] = 1.3
    hdr["SNR"] = 40
    hdr["NVISITS"] = n_spectra

    # Init hdulist
    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU(header=hdr))

    # Init the key HDU's (flux, error, bitmask)
    # names
    units = [
        "Flux (10^-17 erg/s/cm^2/Ang)", "Err (10^-17 erg/s/cm^2/Ang)", "Mask"
    ] + (["test"] * 7)
    for i in range(10):
        hdu = fits.ImageHDU(data=np.random.random((n_spectra, 10)))
        hdu.header["NAXIS"] = 2
        hdu.header["NAXIS1"] = 10
        hdu.header["NAXIS2"] = n_spectra
        hdu.header["CDELT1"] = 6e-06
        hdu.header["CRVAL1"] = 4.179
        hdu.header["BUNIT"] = units[i]
        hdu.name = f"apstar{i}"
        hdulist.append(hdu)

    return hdulist


def apVisit_HDUList():
    """Mock an apVisit HDUList"""
    np.random.seed(20)

    # init primary hdu header
    hdr = fits.Header()
    hdr["FOOBAR"] = "barfoo"
    hdr["SURVEY"] = "SDSS-V"
    hdr["MJD5"] = 99999
    hdr["DATE-OBS"] = "1970-01-01"

    # Init hdulist
    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU(header=hdr))

    # Init the key HDU's (flux, error, bitmask, spectral)
    for i in range(4):
        if i == 3:
            hdu = fits.ImageHDU(data=np.array([
                np.arange(1, 11, 1),
                np.arange(11, 21, 1),
                np.arange(21, 31, 1)
            ]))
            hdu.header["BUNIT"] = "Wavelength (Ang)"

        else:
            hdu = fits.ImageHDU(data=np.random.random((3, 10)))
            hdu.header["BUNIT"] = "Flux (10^-17 erg/s/cm^2/Ang)"
        hdulist.append(hdu)

    return hdulist


def spec_HDUList(n_spectra):
    """Mock an BOSS spec HDUList of n_spectra spectra + 1 coadd."""
    np.random.seed(20)

    # init primary hdu header
    hdr = fits.Header()
    hdr["OBSERVAT"] = "APO"
    hdr["TELESCOP"] = "SDSS 2.5-M"
    hdr["FOOBAR"] = "barfoo"
    hdr["MJD5"] = 99999
    hdr["DATE-OBS"] = "1970-01-01"

    # Init hdulist
    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU(header=hdr))

    # Init the key HDU's (flux, error, bitmask, spectral)
    names = ["COADD", "SPALL", "ZALL", "ZLINE"]
    for i in range(4):
        hdu = fits.BinTableHDU.from_columns([
            fits.Column(name="FLUX", format="E", array=np.random.random(10)),
            fits.Column(name="LOGLAM",
                        format="E",
                        array=np.random.random(10).sort()),
            fits.Column(name="IVAR", format="E", array=np.random.random(10)),
        ])
        hdu.name = names[i]
        hdulist.append(hdu)
    for i in range(n_spectra):
        hdu = fits.BinTableHDU.from_columns([
            fits.Column(name="LOGLAM",
                        format="E",
                        array=np.random.random(10).sort()),
            fits.Column(name="FLUX", format="E", array=np.random.random(10)),
            fits.Column(name="IVAR", format="E", array=np.random.random(10)),
        ])
        hdu.name = f"spectrum{i}"
        hdulist.append(hdu)

    return hdulist


# TEST MWM loaders
@pytest.mark.parametrize(
    "file_obj,hdu, format",
    [
        ("mwmVisit-temp", 3, "visit"),
        ("mwmStar-temp", 4, "star"),
    ],
)
def test_mwm_1d(file_obj, hdu, format):
    """Test mwm Spectrum1D loader"""
    # TODO make work with dummy file
    tmpfile = str(file_obj) + ".fits"
    mwm_HDUList(0, format).writeto(tmpfile, overwrite=True)

    data = Spectrum1D.read(tmpfile, hdu=hdu)
    assert isinstance(data, Spectrum1D)


@pytest.mark.parametrize(
    "file_obj, format",
    [
        ("mwmVisit-temp", "visit"),
        ("mwmStar-temp", "star"),
    ],
)
def test_mwm_list(file_obj, format):
    """Test mwm SpectrumList loader"""
    # TODO make work with dummy file
    tmpfile = str(file_obj) + ".fits"
    mwm_HDUList(0, format).writeto(tmpfile, overwrite=True)

    data = SpectrumList.read(tmpfile, format="SDSS-V mwm multi")
    assert isinstance(data, SpectrumList)
    for i in range(len(data)):
        assert isinstance(data[i], Spectrum1D)


@pytest.mark.parametrize(
    "file_obj,n_spectra",
    [
        ("spec-temp", 1),
        ("spec-temp", 5),
    ],
)
def test_spec_1d(file_obj, n_spectra):
    """Test BOSS spec loader"""
    tmpfile = str(file_obj) + ".fits"
    spec_HDUList(n_spectra).writeto(tmpfile, overwrite=True)

    idxs = [1] + list(np.arange(5, 5 + n_spectra, 1))

    for i in idxs:
        data = Spectrum1D.read(tmpfile, hdu=i)
        assert isinstance(data, Spectrum1D)
        assert len(data.flux.value) == 10
        assert data.flux.unit == Unit("1e-17 erg / (s cm2 Angstrom)")

        assert len(data.spectral_axis.value) == 10
        assert data.spectral_axis.unit == Angstrom

        assert data[i].meta["header"].get("foobar") == "barfoo"


@pytest.mark.parametrize(
    "file_obj,n_spectra",
    [
        ("spec-temp", 0),
        ("spec-temp", 5),
    ],
)
def test_spec_list(file_obj, n_spectra):
    """Test BOSS SpectrumList loader"""
    tmpfile = str(file_obj) + ".fits"
    spec_HDUList(n_spectra).writeto(tmpfile, overwrite=True)

    data = SpectrumList.read(tmpfile, format="SDSS-V spec multi")
    assert isinstance(data, SpectrumList)
    assert len(data) == n_spectra + 1
    for i in range(n_spectra):
        assert len(data[i].flux.value) == 10
        assert data[i].flux.unit == Unit("1e-17 erg / (s cm2 Angstrom)")
        assert len(data[i].spectral_axis.value) == 10
        assert data[i].spectral_axis.unit == Angstrom
        assert data[i].meta["header"].get("foobar") == "barfoo"


@pytest.mark.parametrize(
    "file_obj,hdu",
    [
        ("spec-temp", 2),
        ("spec-temp", 3),
        ("spec-temp", 4),
    ],
    ids=["Fail on HDU2", "Fail on HDU3", "Fail on HDU4"],
)
def test_spec_1d_fail_hdu(file_obj, hdu):
    """Test if fail on reading HDU2, HDU3, or HDU4 (non-spectra)"""
    tmpfile = str(file_obj) + ".fits"
    spec_HDUList(5).writeto(tmpfile, overwrite=True)

    with pytest.raises(ValueError):
        Spectrum1D.read(tmpfile, hdu=hdu)


@pytest.mark.parametrize(
    "file_obj,idx",
    [
        ("apStar-temp", 1),
        ("apStar-temp", 4),
    ],
)
def test_apStar_1D(file_obj, idx):
    tmpfile = str(file_obj) + ".fits"
    apStar_HDUList(6).writeto(tmpfile, overwrite=True)

    data = Spectrum1D.read(tmpfile, idx=idx)
    assert isinstance(data, Spectrum1D)
    assert len(data.flux.value) == 10
    assert data.flux.unit == Unit("1e-17 erg / (s cm2 Angstrom)")

    assert len(data.spectral_axis.value) == 10
    assert data.spectral_axis.unit == Angstrom

    assert data.meta["header"].get("foobar") == "barfoo"


@pytest.mark.parametrize(
    "file_obj,n_spectra",
    [
        ("apStar-temp", 3),
        ("apStar-temp", 6),
    ],
)
def test_apStar_list(file_obj, n_spectra):
    tmpfile = str(file_obj) + ".fits"
    apStar_HDUList(n_spectra).writeto(tmpfile, overwrite=True)

    data = SpectrumList.read(tmpfile, format="SDSS-V apStar multi")
    assert isinstance(data, SpectrumList)
    assert len(data) == n_spectra
    for i in range(len(data)):
        assert isinstance(data[i], Spectrum1D)
        assert len(data[i].flux.value) == 10
        assert data[i].flux.unit == Unit("1e-17 erg / (s cm2 Angstrom)")
        assert len(data[i].spectral_axis.value) == 10
        assert data[i].meta["header"].get("foobar") == "barfoo"


@pytest.mark.parametrize(
    "file_obj",
    [
        ("apStar-temp"),
    ],
)
def test_apStar_fail_list(file_obj):
    """Test if fail on reading 1D apStar as list"""
    tmpfile = str(file_obj) + ".fits"
    apStar_HDUList(1).writeto(tmpfile, overwrite=True)

    with pytest.raises(ValueError):
        SpectrumList.read(tmpfile, format="SDSS-V apStar multi")


@pytest.mark.parametrize(
    "file_obj",
    [
        ("apVisit-temp"),
    ],
)
def test_apVisit_1D(file_obj):
    tmpfile = str(file_obj) + ".fits"
    apVisit_HDUList().writeto(tmpfile, overwrite=True)

    data = Spectrum1D.read(tmpfile)
    assert isinstance(data, Spectrum1D)
    assert np.array_equal(data.spectral_axis.value, np.arange(1, 31, 1))
    assert len(data.flux.value) == 30
    assert data.meta["header"].get("foobar") == "barfoo"


@pytest.mark.parametrize(
    "file_obj",
    [
        ("apVisit-temp"),
    ],
)
def test_apVisit_list(file_obj):
    tmpfile = str(file_obj) + ".fits"
    apVisit_HDUList().writeto(tmpfile, overwrite=True)

    data = SpectrumList.read(tmpfile, format="SDSS-V apVisit multi")
    assert isinstance(data, SpectrumList)
    assert len(data) == 3
