from astropy.nddata import VarianceUncertainty, StdDevUncertainty
import astropy.units as u

from .conftest import remote_access
from .. import SpectrumList
from ..io.default_loaders import dc_common as loaders

REMOTE_ID = "4059032"

GAMA_2QZ_TEST_FILENAME = "J113606.3+001155a.fit"
GAMA_2SLAQ_QSO_TEST_FILENAME = "J091726.21+003424.0_a14_040423.fit"
GAMA_GAMA_LT_TEST_FILENAME = "LTF_09_1128_0555.fit"
GAMA_GAMA_TEST_FILENAME = "G12_Y2_009_044.fit"
GAMA_MGC_TEST_FILENAME = "MGC23320.fit"
GAMA_WIGGLEZ_TEST_FILENAME = "Spectrum-195254.fit"
OZDES_TEST_FILENAME = "OzDES-DR2_04720.fits"
WIGGLEZ_TEST_FILENAME = "wig206635.fits"

OZDES_CONFIG = {
    "hdus": {
        "0": {"purpose": "combined_science"},
        "1": {"purpose": "combined_error_variance"},
        "2": {"purpose": "skip"},
        "cycle": {
            "0": {"purpose": "science"},
            "1": {"purpose": "error_variance"},
            "2": {"purpose": "skip"},
        },
    },
    "units": None,
    "wcs": None,
    "all_standard_units": True,
    "all_keywords": False,
    "valid_wcs": True,
}
GAMA_2QZ_CONFIG = {
    "hdus": None,
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CD1_1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": True,
    "valid_wcs": False,
}
GAMA_2SLAQ_QSO_CONFIG = {
    "hdus": None,
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CDELT1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": True,
    "valid_wcs": False,
}
GAMA_LT_CONFIG = {
    "hdus": {"0": {"purpose": "science"}, },
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX",
        "pixel_reference_point_value_keyword": "CRVAL",
        "pixel_width_keyword": "CDELT",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": False,
    "valid_wcs": False,
}
GAMA_WIGGLEZ_CONFIG = {
    "hdus": {
        "0": {"purpose": "science"},
        "1": {"purpose": "error_variance"},
        "2": {"purpose": "skip"},
    },
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CDELT1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": False,
    "valid_wcs": False,
}
GAMA_GAMA_CONFIG = {
    "hdu": {
        "1": {
            "purpose": "science",
            "units": {"flux_unit": "10^-17 erg/s/cm^2/A"},
        },
        "2": {"purpose": "error_stdev"},
        "3": {"purpose": "unreduced_science"},
        "4": {"purpose": "unreduced_error_stdev"},
        "5": {"purpose": "sky"},
    },
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CD1_1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": False,
    "valid_wcs": False,
}
GAMA_MGC_CONFIG = {
    "hdu": None,
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CD1_1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": True,
    "valid_wcs": False,
}


class TestSingleSplit:
    @remote_access([{'id': "4460981", 'filename': WIGGLEZ_TEST_FILENAME}])
    def test_wigglez(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path, format="WiggleZ")

        assert len(spectra) == 1
        assert spectra[0].flux.unit == u.Unit("1e-16 erg / (A cm2 s)")
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)
        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @remote_access([{'id': "4460981", 'filename': WIGGLEZ_TEST_FILENAME}])
    def test_wigglez_guess(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path)

        assert len(spectra) == 1
        assert spectra[0].flux.unit == u.Unit("1e-16 erg / (A cm2 s)")
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)
        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': OZDES_TEST_FILENAME}])
    def test_ozdes(self, remote_data_path):
        spectra = SpectrumList.read(
            remote_data_path, format=loaders.SINGLE_SPLIT_LABEL, **OZDES_CONFIG
        )

        # The test file has the combined obs, and 4 other sets
        assert len(spectra) == 5

        assert spectra[0].flux.unit == u.count / u.Angstrom
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)
        assert spectra[0].meta.get("label") is None
        assert spectra[0].meta.get("header") is not None

        assert spectra[1].flux.unit == u.count / u.Angstrom
        assert spectra[1].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[1].uncertainty, VarianceUncertainty)
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None

        assert spectra[2].flux.unit == u.count / u.Angstrom
        assert spectra[2].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[2].uncertainty, VarianceUncertainty)
        assert spectra[2].meta.get("label") is not None
        assert spectra[2].meta.get("header") is not None

        assert spectra[3].flux.unit == u.count / u.Angstrom
        assert spectra[3].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[3].uncertainty, VarianceUncertainty)
        assert spectra[3].meta.get("label") is not None
        assert spectra[3].meta.get("header") is not None

        assert spectra[4].flux.unit == u.count / u.Angstrom
        assert spectra[4].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[4].uncertainty, VarianceUncertainty)
        assert spectra[4].meta.get("label") is not None
        assert spectra[4].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': OZDES_TEST_FILENAME}])
    def test_ozdes_named_loader(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path, format="OzDES")

        # The test file has the combined obs, and 4 other sets
        assert len(spectra) == 5

        assert spectra[0].flux.unit == u.count / u.Angstrom
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)
        assert spectra[0].meta.get("label") is None
        assert spectra[0].meta.get("header") is not None

        assert spectra[1].flux.unit == u.count / u.Angstrom
        assert spectra[1].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[1].uncertainty, VarianceUncertainty)
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None

        assert spectra[2].flux.unit == u.count / u.Angstrom
        assert spectra[2].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[2].uncertainty, VarianceUncertainty)
        assert spectra[2].meta.get("label") is not None
        assert spectra[2].meta.get("header") is not None

        assert spectra[3].flux.unit == u.count / u.Angstrom
        assert spectra[3].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[3].uncertainty, VarianceUncertainty)
        assert spectra[3].meta.get("label") is not None
        assert spectra[3].meta.get("header") is not None

        assert spectra[4].flux.unit == u.count / u.Angstrom
        assert spectra[4].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[4].uncertainty, VarianceUncertainty)
        assert spectra[4].meta.get("label") is not None
        assert spectra[4].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_2QZ_TEST_FILENAME}])
    def test_gama_2qz(self, remote_data_path):
        spectra = SpectrumList.read(
            remote_data_path, format=loaders.SINGLE_SPLIT_LABEL,
            **GAMA_2QZ_CONFIG
        )
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_2QZ_TEST_FILENAME}])
    def test_gama_2qz_named_loader(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path, format="GAMA-2QZ")
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_2QZ_TEST_FILENAME}])
    def test_gama_2qz_guess(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path)
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_2SLAQ_QSO_TEST_FILENAME}])
    def test_2slaq_qso(self, remote_data_path):
        spectra = SpectrumList.read(
            remote_data_path, format=loaders.SINGLE_SPLIT_LABEL,
            **GAMA_2SLAQ_QSO_CONFIG
        )
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_2SLAQ_QSO_TEST_FILENAME}])
    def test_2slaq_qso_named_loader(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path, format="GAMA-2SLAQ-QSO")
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_2SLAQ_QSO_TEST_FILENAME}])
    def test_2slaq_qso_normalised(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path)
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_GAMA_LT_TEST_FILENAME}])
    def test_gama_lt(self, remote_data_path):
        spectra = SpectrumList.read(
            remote_data_path,
            format=loaders.SINGLE_SPLIT_LABEL,
            **GAMA_LT_CONFIG
        )
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert spectra[0].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_GAMA_LT_TEST_FILENAME}])
    def test_gama_lt_named_loader(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path, format="GAMA-LT")
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert spectra[0].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_GAMA_LT_TEST_FILENAME}])
    def test_gama_lt_guess(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path)
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert spectra[0].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_WIGGLEZ_TEST_FILENAME}])
    def test_gama_wigglez(self, remote_data_path):
        spectra = SpectrumList.read(
            remote_data_path, format=loaders.SINGLE_SPLIT_LABEL,
            **GAMA_WIGGLEZ_CONFIG
        )
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_WIGGLEZ_TEST_FILENAME}])
    def test_gama_wigglez_named_loader(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path, format="GAMA-WiggleZ")
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_WIGGLEZ_TEST_FILENAME}])
    def test_gama_wigglez_guess(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path)
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None


class TestMultilineSingle:
    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_GAMA_TEST_FILENAME}])
    def test_gama_gama(self, remote_data_path):
        spectra = SpectrumList.read(
            remote_data_path, format=loaders.MULTILINE_SINGLE_LABEL,
            **GAMA_GAMA_CONFIG
        )
        assert len(spectra) == 3

        assert spectra[0].flux.unit == u.Unit("10^-17 erg/s/cm^2/A")
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[1].flux.unit == u.count
        assert spectra[1].spectral_axis.unit == u.Angstrom
        assert spectra[2].flux.unit == u.count
        assert spectra[2].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, StdDevUncertainty)
        assert isinstance(spectra[1].uncertainty, StdDevUncertainty)
        assert spectra[2].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None
        assert spectra[2].meta.get("label") is not None
        assert spectra[2].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_GAMA_TEST_FILENAME}])
    def test_gama_gama_named_loader(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path, format="GAMA")
        assert len(spectra) == 3

        assert spectra[0].flux.unit == u.Unit("10^-17 erg/s/cm^2/A")
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[1].flux.unit == u.count
        assert spectra[1].spectral_axis.unit == u.Angstrom
        assert spectra[2].flux.unit == u.count
        assert spectra[2].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, StdDevUncertainty)
        assert isinstance(spectra[1].uncertainty, StdDevUncertainty)
        assert spectra[2].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None
        assert spectra[2].meta.get("label") is not None
        assert spectra[2].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_GAMA_TEST_FILENAME}])
    def test_gama_gama_guess(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path)
        assert len(spectra) == 3

        assert spectra[0].flux.unit == u.Unit("10^-17 erg/s/cm^2/A")
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[1].flux.unit == u.count
        assert spectra[1].spectral_axis.unit == u.Angstrom
        assert spectra[2].flux.unit == u.count
        assert spectra[2].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, StdDevUncertainty)
        assert isinstance(spectra[1].uncertainty, StdDevUncertainty)
        assert spectra[2].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None
        assert spectra[2].meta.get("label") is not None
        assert spectra[2].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_MGC_TEST_FILENAME}])
    def test_gama_mgc(self, remote_data_path):
        spectra = SpectrumList.read(
            remote_data_path, format=loaders.MULTILINE_SINGLE_LABEL,
            **GAMA_MGC_CONFIG
        )
        assert len(spectra) == 2

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[1].flux.unit == u.count
        assert spectra[1].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, StdDevUncertainty)
        assert spectra[1].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_MGC_TEST_FILENAME}])
    def test_gama_mgc_named_loader(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path, format="GAMA-MGC")
        assert len(spectra) == 2

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[1].flux.unit == u.count
        assert spectra[1].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, StdDevUncertainty)
        assert spectra[1].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None

    @remote_access([{'id': REMOTE_ID, 'filename': GAMA_MGC_TEST_FILENAME}])
    def test_gama_mgc_guess(self, remote_data_path):
        spectra = SpectrumList.read(remote_data_path)
        assert len(spectra) == 2

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[1].flux.unit == u.count
        assert spectra[1].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, StdDevUncertainty)
        assert spectra[1].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None
