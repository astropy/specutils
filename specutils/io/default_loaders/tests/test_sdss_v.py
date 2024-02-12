import numpy as np
import pytest

from astropy.io import fits
from astropy.units import Unit, Angstrom

from specutils import Spectrum1D, SpectrumList


def generate_apogee_hdu(observatory="APO", with_wl=True, datasum="0"):
    wl = (10**(4.179 + 6e-6 * np.arange(8575))).reshape((1, -1))
    flux = np.zeros_like(wl)
    ivar = np.zeros_like(wl)
    pixel_flags = np.zeros_like(wl)
    continuum = np.zeros_like(wl)
    nmf_rectified_model_flux = np.zeros_like(wl)

    columns = [
        fits.Column(name="spectrum_pk_id", array=[159783564], format="K"),
        fits.Column(name="release", array=[b"sdss5"], format="5A"),
        fits.Column(name="filetype", array=[b"apStar"], format="6A"),
        fits.Column(name="v_astra", array=[b"0.5.0"], format="5A"),
        fits.Column(name="healpix", array=[3], format="J"),
        fits.Column(name="sdss_id", array=[42], format="K"),
        fits.Column(name="apred", array=[b"1.2"], format="3A"),
        fits.Column(name="obj", array=[b"2M19534321+6705175"], format="18A"),
        fits.Column(name="telescope", array=[b"apo25m"], format="6A"),
        fits.Column(name="min_mjd", array=[59804], format="J"),
        fits.Column(name="max_mjd", array=[59866], format="J"),
        fits.Column(name="n_entries", array=[-1], format="J"),
        fits.Column(name="n_visits", array=[5], format="J"),
        fits.Column(name="n_good_visits", array=[5], format="J"),
        fits.Column(name="n_good_rvs", array=[5], format="J"),
        fits.Column(name="snr", array=[46.56802], format="E"),
        fits.Column(name="mean_fiber", array=[256.0], format="E"),
        fits.Column(name="std_fiber", array=[0.0], format="E"),
        fits.Column(name="spectrum_flags", array=[1048576], format="J"),
        fits.Column(name="v_rad", array=[-56.7284381], format="E"),
        fits.Column(name="e_v_rad", array=[5.35407624], format="E"),
        fits.Column(name="std_v_rad", array=[10.79173857], format="E"),
        fits.Column(name="median_e_v_rad", array=[16.19418386], format="E"),
        fits.Column(name="doppler_teff", array=[7169.0107], format="E"),
        fits.Column(name="doppler_e_teff", array=[9.405238], format="E"),
        fits.Column(name="doppler_logg", array=[2.981389], format="E"),
        fits.Column(name="doppler_e_logg", array=[0.01916536], format="E"),
        fits.Column(name="doppler_fe_h", array=[-1.20532212], format="E"),
        fits.Column(name="doppler_e_fe_h", array=[0.0093738], format="E"),
        fits.Column(name="doppler_rchi2", array=[1.1424173], format="E"),
        fits.Column(name="doppler_flags", array=[0], format="J"),
        fits.Column(name="xcorr_v_rad", array=[np.nan], format="E"),
        fits.Column(name="xcorr_v_rel", array=[np.nan], format="E"),
        fits.Column(name="xcorr_e_v_rel", array=[np.nan], format="E"),
        fits.Column(name="ccfwhm", array=[np.nan], format="E"),
        fits.Column(name="autofwhm", array=[np.nan], format="E"),
        fits.Column(name="n_components", array=[1], format="J"),
    ]
    if with_wl:
        columns.append(
            fits.Column(name="wavelength",
                        array=wl,
                        format="8575E",
                        dim="(8575)"))
    columns += [
        fits.Column(name="flux", array=flux, format="8575E", dim="(8575)"),
        fits.Column(name="ivar", array=ivar, format="8575E", dim="(8575)"),
        fits.Column(name="pixel_flags",
                    array=pixel_flags,
                    format="8575E",
                    dim="(8575)"),
        fits.Column(name="continuum",
                    array=continuum,
                    format="8575E",
                    dim="(8575)"),
        fits.Column(
            name="nmf_rectified_model_flux",
            array=nmf_rectified_model_flux,
            format="8575E",
            dim="(8575)",
        ),
        fits.Column(name="nmf_rchi2", array=[2.3391197], format="E"),
        fits.Column(name="nmf_flags", array=[0], format="J"),
    ]
    header = fits.Header(cards=[
        ("EXTNAME", f"APOGEE/{observatory}", ""),
        ("OBSRVTRY", observatory, None),
        ("INSTRMNT", "APOGEE", None),
        ("CRVAL", 4.179, None),
        ("CDELT", 6e-6, None),
        ("CTYPE", "LOG-LINEAR", None),
        ("CUNIT", "Angstrom (Vacuum)"),
        ("CRPIX", 1, None),
        ("DC-FAG", 1, None),
        ("NPIXELS", 8575, None),
        ("DATASUM", datasum, "data unit checksum updated 2023-11-13T03:21:47"),
    ])

    return fits.BinTableHDU.from_columns(columns, header=header)


def generate_boss_hdu(observatory="APO", with_wl=True, datasum="0"):
    wl = (10**(3.5523 + 1e-4 * np.arange(4648))).reshape((1, -1))
    flux = ivar = continuum = pixel_flags = nmf_rectified_model_flux = np.zeros_like(
        wl)
    columns = [
        fits.Column(name="spectrum_pk_id", array=[0], format="K"),
        fits.Column(name="release", array=["sdss5"], format="5A"),
        fits.Column(name="filetype", array=["specFull"], format="7A"),
        fits.Column(name="v_astra", array=["0.5.0"], format="5A"),
        fits.Column(name="healpix", array=[34], format="J"),
        fits.Column(name="sdss_id", array=[42], format="K"),
        fits.Column(name="run2d", array=["6_1_2"], format="6A"),
        fits.Column(name="telescope", array=["apo25m"], format="6A"),
        fits.Column(name="min_mjd", array=[54], format="J"),
        fits.Column(name="max_mjd", array=[488], format="J"),
        fits.Column(name="n_visits", array=[1], format="J"),
        fits.Column(name="n_good_visits", array=[1], format="J"),
        fits.Column(name="n_good_rvs", array=[1], format="J"),
        fits.Column(name="v_rad", array=[0], format="E"),
        fits.Column(name="e_v_rad", array=[1], format="E"),
        fits.Column(name="std_v_rad", array=[1], format="E"),
        fits.Column(name="median_e_v_rad", array=[3], format="E"),
        fits.Column(name="xcsao_teff", array=[5000], format="E"),
        fits.Column(name="xcsao_e_teff", array=[10], format="E"),
        fits.Column(name="xcsao_logg", array=[4], format="E"),
        fits.Column(name="xcsao_e_logg", array=[3], format="E"),
        fits.Column(name="xcsao_fe_h", array=[0], format="E"),
        fits.Column(name="xcsao_e_fe_h", array=[5], format="E"),
        fits.Column(name="xcsao_meanrxc", array=[0], format="E"),
        fits.Column(name="snr", array=[50], format="E"),
        fits.Column(name="gri_gaia_transform_flags", array=[1], format="J"),
        fits.Column(name="zwarning_flags", array=[0], format="J"),
    ]
    if with_wl:
        columns.append(
            fits.Column(name="wavelength",
                        array=wl,
                        format="4648E",
                        dim="(4648)"))
    columns += [
        fits.Column(name="flux", array=flux, format="4648E", dim="(4648)"),
        fits.Column(name="ivar", array=ivar, format="4648E", dim="(4648)"),
        fits.Column(name="pixel_flags",
                    array=pixel_flags,
                    format="4648E",
                    dim="(4648)"),
        fits.Column(name="continuum",
                    array=continuum,
                    format="4648E",
                    dim="(4648)"),
        fits.Column(
            name="nmf_rectified_model_flux",
            array=nmf_rectified_model_flux,
            format="4648E",
            dim="(4648)",
        ),
        fits.Column(name="nmf_rchi2", array=[5], format="E"),
        fits.Column(name="nmf_flags", array=[0], format="J"),
    ]
    header = fits.Header(cards=[
        ("EXTNAME", f"BOSS/{observatory}", ""),
        ("OBSRVTRY", observatory, None),
        ("INSTRMNT", "BOSS", None),
        ("CRVAL", 3.5523, None),
        ("CDELT", 1e-4, None),
        ("CTYPE", "LOG-LINEAR", None),
        ("CUNIT", "Angstrom (Vacuum)"),
        ("CRPIX", 1, None),
        ("DC-FAG", 1, None),
        ("NPIXELS", 4648, None),
        ("DATASUM", datasum, "data unit checksum updated 2023-11-13T03:21:47"),
    ])

    return fits.BinTableHDU.from_columns(columns, header=header)


def fake_primary_hdu():
    return fits.PrimaryHDU(header=fits.Header(cards=[
        ("SIMPLE", True, "conforms to FITS standard"),
        ("BITPIX", 8, "array data type"),
        ("NAXIS", 0, "number of array dimensions"),
        ("EXTEND", True, ""),
        ("", "", ""),
        ("", "Metadata", ""),
        ("", "", ""),
        ("V_ASTRA", "0.5.0", "Astra version"),
        (
            "CREATED",
            "23-11-13 10:21:47",
            "File creation time (UTC %y-%m-%d %H:%M:%S)",
        ),
        ("", "", ""),
        ("", "Identifiers", ""),
        ("", "", ""),
        ("SDSS_ID", 34, "SDSS-5 unique identifier"),
        ("APOGEEID", "", "SDSS-4 DR17 APOGEE identifier"),
        ("GAIA2_ID", 23423, "Gaia DR2 source identifier"),
        ("GAIA3_ID", 23423, "Gaia DR3 source identifier"),
        ("TIC_ID", 453, "TESS Input Catalog (v8) identifier"),
        ("HEALPIX", 7769, "HEALPix (128 side)"),
        ("", "", ""),
        ("", "Targeting Provenance", ""),
        ("", "", ""),
        ("CARTON_0", "", "Highest priority carton name"),
        ("LEAD", "tic_v8", "Lead catalog used for cross-match"),
        ("VER_ID", 25, "SDSS catalog version for targeting"),
        (
            "CAT_ID",
            27021597842679140,
            "Catalog identifier used to target the source",
        ),
        ("CAT_ID21", 4283338543, "Catalog identifier (v21; v0.0)"),
        ("CAT_ID25", 27021597842679140, "Catalog identifier (v25; v0.5)"),
        ("CAT_ID31", 63050395058384891, "Catalog identifier (v31; v1.0)"),
        ("N_ASSOC", 1, "SDSS_IDs associated with this CATALOGID"),
        ("N_NEIGH", 0, 'Sources within 3" and G_MAG < G_MAG_source + 5'),
        ("", "", ""),
        ("", "Astrometry", ""),
        ("", "", ""),
        ("RA", 198.4346, "Right ascension [deg]"),
        ("DEC", 67.110405, "Declination [deg]"),
        ("L", 39.57304281336641, "Galactic longitude [deg]"),
        ("B", 69.09476252445475, "Galactic latitude [deg]"),
        ("PLX", 6.2856047, "Parallax [mas]"),
        ("E_PLX", 1.04624034, "Error on parallax [mas]"),
        ("PMRA", -3.6440573, "Proper motion in RA [mas/yr]"),
        ("E_PMRA", 6.05078333, "Error on proper motion in RA [mas/yr]"),
        ("PMDE", -8.818821, "Proper motion in DEC [mas/yr]"),
        ("E_PMDE", 1.05103116, "Error on proper motion in DEC [mas/yr]"),
        ("V_RAD", "NaN", "Gaia radial velocity [km/s]"),
        ("E_V_RAD", "NaN", "Error on Gaia radial velocity [km/s]"),
        ("", "", ""),
        ("", "Gaia Photometry", ""),
        ("", "", ""),
        ("G_MAG", 5.212288, "Gaia DR3 mean G band magnitude [mag]"),
        ("BP_MAG", 1.417452, "Gaia DR3 mean BP band magnitude [mag]"),
        ("RP_MAG", 13.138654, "Gaia DR3 mean RP band magnitude [mag]"),
        ("", "", ""),
        ("", "2MASS Photometry", ""),
        ("", "", ""),
        ("J_MAG", 3.709, "2MASS J band magnitude [mag]"),
        ("E_J_MAG", 0.026, "Error on 2MASS J band magnitude [mag]"),
        ("H_MAG", 2.893, "2MASS H band magnitude [mag]"),
        ("E_H_MAG", 0.032, "Error on 2MASS H band magnitude [mag]"),
        ("K_MAG", 1.776, "2MASS K band magnitude [mag]"),
        ("E_K_MAG", 0.026, "Error on 2MASS K band magnitude [mag]"),
        ("PH_QUAL", "", "2MASS photometric quality flag"),
        ("BL_FLG", "", "Number of components fit per band (JHK)"),
        ("CC_FLG", "", "Contamination and confusion flag"),
        (
            "COMMENT",
            "See https://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec2_2a.html",
            "",
        ),
        ("", "", ""),
        ("", "unWISE Photometry", ""),
        ("", "", ""),
        ("W1_MAG", 1.615836388268159, "W1 magnitude"),
        ("E_W1_MAG", 3.001952228308748062, "Error on W1 magnitude"),
        ("W1_FLUX", 954.997, "W1 flux [Vega nMgy]"),
        ("W1_DFLUX", 1.1017, "Error on W1 flux [Vega nMgy]"),
        ("W1_FRAC", 0.976102, "Fraction of W1 flux from this object"),
        ("W2_MAG", 12.20587252800215, "W2 magnitude [Vega]"),
        ("E_W2_MAG", 0.04519932356304733, "Error on W2 magnitude"),
        ("W2_FLUX", 868.906, "W2 flux [Vega nMgy]"),
        ("W2_DFLUX", 3.172016, "Error on W2 flux [Vega nMgy]"),
        ("W2_FRAC", 0.9741495, "Fraction of W2 flux from this object"),
        ("W1UFLAGS", 0, "unWISE flags for W1"),
        ("W2UFLAGS", 0, "unWISE flags for W2"),
        ("W1AFLAGS", 0, "Additional flags for W1"),
        ("W2AFLAGS", 0, "Additional flags for W2"),
        ("COMMENT", "See https://catalog.unwise.me/catalogs.html", ""),
        ("", "", ""),
        ("", "GLIMPSE Photometry", ""),
        ("", "", ""),
        ("MAG4_5", "", "IRAC band 4.5 micron magnitude [mag]"),
        ("D4_5M", "", "Error on IRAC band 4.5 micron magnitude [mag]"),
        ("RMS_F4_5", "", "RMS deviations from final flux [mJy]"),
        ("SQF_4_5", 0, "Source quality flag for IRAC band 4.5 micron"),
        ("MF4_5", 0, "Flux calculation method flag"),
        ("CSF", 0, "Close source flag"),
        (
            "COMMENT",
            "See https://irsa.ipac.caltech.edu/data/SPITZER/GLIMPSE/gator_docs/",
            "",
        ),
        ("", "", ""),
        ("", "Observations Summary", ""),
        ("", "", ""),
        ("N_BOSS", 1, "Number of BOSS visits"),
        ("B_MINMJD", 59804, "Minimum MJD of BOSS visits"),
        ("N_MAXMJD", 59866, "Maximum MJD of BOSS visits"),
        ("N_APOGEE", 1, "Number of APOGEE visits"),
        ("A_MINMJD", 59804, "Minimum MJD of APOGEE visits"),
        ("A_MAXMJD", 59866, "Maximum MJD of APOGEE visits"),
        ("", "", ""),
        ("", "Data Integrity", ""),
        ("", "", ""),
        (
            "CHECKSUM",
            "UKTRaHQOWHQOaHQO",
            "HDU checksum updated 2023-11-13T03:21:47",
        ),
        ("DATASUM", "0", "data unit checksum updated 2023-11-13T03:21:47"),
        ("", "", ""),
        ("", "HDU Descriptions", ""),
        ("", "", ""),
        ("COMMENT", "HDU 0: Summary information only", ""),
        ("COMMENT", "HDU 1: BOSS spectra from Apache Point Observatory", ""),
        ("COMMENT", "HDU 2: BOSS spectra from Las Campanas Observatory", ""),
        ("COMMENT", "HDU 3: APOGEE spectra from Apache Point Observatory", ""),
        ("COMMENT", "HDU 4: APOGEE spectra from Las Campanas Observatory", ""),
    ]))


def mwm_HDUList(hduflags, with_wl):
    hdulist = [fake_primary_hdu()]
    for i, flag in enumerate(hduflags):
        obs = ["APO", "LCO"]
        if i <= 1:
            hdulist.append(
                generate_boss_hdu(obs[i % 2],
                                  with_wl=with_wl,
                                  datasum=str(flag)))
        else:
            hdulist.append(
                generate_apogee_hdu(obs[i % 2],
                                    with_wl=with_wl,
                                    datasum=str(flag)))

    print(hdulist)
    return fits.HDUList(hdulist)


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
            fits.Column(name="AND_MASK",
                        format="E",
                        array=np.random.random(10)),
            fits.Column(name="OR_MASK", format="E",
                        array=np.random.random(10)),
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
            fits.Column(name="AND_MASK",
                        format="E",
                        array=np.random.random(10)),
            fits.Column(name="OR_MASK", format="E",
                        array=np.random.random(10)),
        ])
        hdu.name = f"spectrum{i}"
        hdulist.append(hdu)

    return hdulist


# TEST MWM loaders
@pytest.mark.parametrize(
    "file_obj, hdu, with_wl, hduflags",
    [
        ("mwm-temp", None, False, [0, 0, 1, 0]),
        ("mwm-temp", 3, False, [0, 0, 1, 0]),
        ("mwm-temp", None, True, [0, 1, 1, 0]),
        ("mwm-temp", 2, True, [0, 1, 1, 0]),
    ],
)
def test_mwm_1d(file_obj, hdu, with_wl, hduflags):
    """Test mwm Spectrum1D loader"""
    tmpfile = str(file_obj) + ".fits"
    mwm_HDUList(hduflags, with_wl).writeto(tmpfile, overwrite=True)

    if hdu is None:
        data = Spectrum1D.read(tmpfile)
    else:
        data = Spectrum1D.read(tmpfile, hdu=hdu)
    assert isinstance(data, Spectrum1D)
    assert isinstance(data.meta["header"], fits.Header)
    if data.meta["instrument"].lower() == "apogee":
        length = 8575
    elif data.meta["instrument"].lower() == "boss":
        length = 4648
    else:
        raise ValueError(
            "INSTRMNT tag in test HDU header is not set properly.")
    assert len(data.spectral_axis.value) == length
    assert len(data.flux.value) == length
    assert data.spectral_axis.unit == Angstrom
    assert data.flux.unit == Unit("1e-17 erg / (s cm2 Angstrom)")


@pytest.mark.parametrize(
    "file_obj, with_wl, hduflags",
    [
        ("mwm-temp", False, [1, 1, 1, 1]),
        ("mwm-temp", True, [0, 1, 0, 1]),
    ],
)
def test_mwm_list(file_obj, with_wl, hduflags):
    """Test mwm SpectrumList loader"""
    tmpfile = str(file_obj) + ".fits"
    mwm_HDUList(hduflags, with_wl).writeto(tmpfile, overwrite=True)

    data = SpectrumList.read(tmpfile, format="SDSS-V mwm multi")
    assert isinstance(data, SpectrumList)
    for i in range(len(data)):
        assert isinstance(data[i], Spectrum1D)
        assert isinstance(data[i].meta["header"], fits.Header)
        if data[i].meta["instrument"].lower() == "apogee":
            length = 8575
        elif data[i].meta["instrument"].lower() == "boss":
            length = 4648
        else:
            raise ValueError(
                "INSTRMNT tag in test HDU header is not set properly.")
        assert len(data[i].spectral_axis.value) == length
        assert len(data[i].flux.value) == length
        assert data[i].spectral_axis.unit == Angstrom
        assert data[i].flux.unit == Unit("1e-17 erg / (s cm2 Angstrom)")


@pytest.mark.parametrize(
    "file_obj, hdu, hduflags",
    [
        ("mwm-temp", 1, [0, 1, 0, 1]),
        ("mwm-temp", 2, [1, 0, 1, 1]),
        ("mwm-temp", 3, [0, 0, 0, 1]),
        ("mwm-temp", 4, [0, 1, 1, 0]),
    ],
)
def test_mwm_1d_fail_spec(file_obj, hdu, hduflags):
    """Test mwm Spectrum1D loader fail on bad spec"""
    tmpfile = str(file_obj) + ".fits"
    mwm_HDUList(hduflags, True).writeto(tmpfile, overwrite=True)
    with pytest.raises(IndexError):
        Spectrum1D.read(tmpfile, hdu=hdu)


@pytest.mark.parametrize(
    "file_obj, with_wl",
    [
        ("mwm-temp", False),
        ("mwm-temp", True),
    ],
)
def test_mwm_1d_fail(file_obj, with_wl):
    """Test mwm Spectrum1D loader fail on empty"""
    tmpfile = str(file_obj) + ".fits"
    mwm_HDUList([0, 0, 0, 0], with_wl).writeto(tmpfile, overwrite=True)

    with pytest.raises(ValueError):
        Spectrum1D.read(tmpfile)


@pytest.mark.parametrize(
    "file_obj, with_wl",
    [
        ("mwm-temp", False),
        ("mwm-temp", True),
    ],
)
def test_mwm_list_fail(file_obj, with_wl):
    """Test mwm SpectrumList loader fail on empty"""
    tmpfile = str(file_obj) + ".fits"
    mwm_HDUList([0, 0, 0, 0], with_wl).writeto(tmpfile, overwrite=True)

    with pytest.raises(ValueError):
        SpectrumList.read(tmpfile, format="SDSS-V mwm multi")


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
        assert len(data.mask) == 10

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
        assert len(data[i].mask) == 10
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
