import astropy.units as u
import numpy as np
import pytest
from astropy.utils.exceptions import AstropyUserWarning

from specutils import Spectrum, SpectrumList
from specutils.io.default_loaders import roman as roman_loaders  # noqa: F401

asdf = pytest.importorskip("asdf")


def _roman_tree(mode='single', optel="GRISM", unit="W m**(-2) nm**(-1)", nsrc=1):
    """Function to create a fake asdf tree for roman spectral data"""
    wl = np.linspace(1.0, 2.0, 5, dtype=float)
    flux = np.arange(5, dtype=float) + 1.0
    ferr = np.full_like(flux, 0.1)

    meta = {
                "optical_element": optel,
                "title": "my_1d_spectra_product",
                "unit_wl": "nm",
                "unit_flux": unit,
            }

    if mode == 'single':
        data = {
                "wl": wl,
                "flux": flux,
                "flux_error": ferr,
            }
    elif mode == 'multi':
        data = {}
        for i in range(nsrc):
            data[str(i)] = {
                "wl": wl,
                "flux": (np.arange(5, dtype=float) + 1.0) * (i + 1),
                "flux_error": np.full_like(flux, 0.1 * (i + 1)),
            }

    return {
        "roman": {
            "meta": meta,
            "data": data
        }
    }

@pytest.fixture()
def roman_single(tmp_path):
    """Fixture to create a single source roman spectral asdf file"""
    path = tmp_path / "roman_test_single.asdf"
    af = asdf.AsdfFile(_roman_tree(mode='single'))
    af.write_to(path)
    yield str(path)
    af.close()


@pytest.fixture()
def roman_multi(tmp_path):
    """Fixture to create a roman ssc 1d_combined spectral file"""
    path = tmp_path / "roman_test_multi.asdf"
    af = asdf.AsdfFile(_roman_tree(mode='multi', nsrc=3))
    af.write_to(path)
    yield str(path)
    af.close()


@pytest.fixture()
def roman_file(tmp_path):
    """Fixture factory to create variations on roman spectral products"""
    def _roman_file(mode='single', optel='GRISM', unit="W m**(-2) nm**(-1)", nsrc=1):
        path = tmp_path / "roman_test.asdf"
        af = asdf.AsdfFile(_roman_tree(mode=mode, optel=optel, unit=unit, nsrc=nsrc))
        af.write_to(path)
        af.close()
        return str(path)
    yield _roman_file


def test_identify_roman_single(roman_single):
    """test we can identify a roman single source spectrum"""
    assert roman_loaders.identify_1d_single_source(None, roman_single) is True
    assert roman_loaders.identify_1d_combined(None, roman_single) is False
    assert roman_loaders.identify_1d_individual(None, roman_single) is False

def test_identify_roman_combined(roman_multi):
    """test we can identify a roman 1d combined spectra file"""
    assert roman_loaders.identify_1d_single_source(None, roman_multi) is False
    assert roman_loaders.identify_1d_combined(None, roman_multi) is True
    assert roman_loaders.identify_1d_individual(None, roman_multi) is False

def test_identify_roman_individual(roman_file):
    """test we can identify a roman 1d individual spectra file"""
    roman_multi = roman_file(mode='multi', unit="DN/s")
    assert roman_loaders.identify_1d_single_source(None, roman_multi) is False
    assert roman_loaders.identify_1d_combined(None, roman_multi) is False
    assert roman_loaders.identify_1d_individual(None, roman_multi) is True

def test_identify_not_spectra(roman_file):
    """test that another roman file is not mis-identified"""
    roman_wfi = roman_file(mode='single', optel='WFI', unit="DN/s")
    assert roman_loaders.identify_1d_single_source(None, roman_wfi) is False
    assert roman_loaders.identify_1d_combined(None, roman_wfi) is False
    assert roman_loaders.identify_1d_individual(None, roman_wfi) is False


def test_roman_1d_spectrum(roman_single):
    """test loading a roman 1d single source spectrum"""
    spec = Spectrum.read(roman_single, format="Roman 1d spectra")

    assert isinstance(spec, Spectrum)
    assert spec.spectral_axis.unit == u.nm
    assert spec.unit == u.Unit("W m**(-2) nm**(-1)")
    assert spec.shape == (5,)
    assert spec.uncertainty is not None
    assert np.allclose(spec.uncertainty.array, 0.1)


def test_roman_1d_combined_load_specific_source(roman_multi):
    """test loading a roman 1d combined spectrum for a specific source"""
    spec = Spectrum.read(roman_multi, format="Roman 1d combined", source="2")
    assert isinstance(spec, Spectrum)
    assert spec.meta["source_id"] == "2"
    assert spec.spectral_axis.unit == u.nm
    assert spec.unit == u.Unit("W m**(-2) nm**(-1)")
    assert spec.shape == (5,)


def test_roman_1d_combined_load_first_source(roman_multi):
    """test loading without a source throws the warning and loads first source"""
    with pytest.warns(AstropyUserWarning, match="source not specified"):
        spec = Spectrum.read(roman_multi, format="Roman 1d combined")

    assert spec.meta["source_id"] == "0"


def test_roman_1d_combined_invalid_source(roman_multi):
    """test loading with a bad source id raises error"""
    with pytest.raises(ValueError, match="Invalid source"):
        Spectrum.read(roman_multi, format="Roman 1d combined", source="nope")


def test_roman_1d_combined_list(roman_multi):
    """test we can load a 1d combined list"""
    speclist = SpectrumList.read(roman_multi, format="Roman 1d combined")

    assert isinstance(speclist, SpectrumList)
    assert len(speclist) == 3
    for i, sp in enumerate(speclist):
        assert isinstance(sp, Spectrum)
        assert sp.meta["source_id"] == str(i)
        assert "source_map" in sp.meta
        assert sp.meta["source_map"][str(i)] == i


def test_roman_1d_individual_list(roman_file):
    """test we can load a 1d individual list"""
    roman_indiv = roman_file(mode='multi', unit="DN/s", nsrc=3)
    speclist = SpectrumList.read(roman_indiv, format="Roman 1d individual")

    assert isinstance(speclist, SpectrumList)
    assert len(speclist) == 3
    assert speclist[0].unit == u.Unit("DN/s")
    assert speclist[0].spectral_axis.unit == u.nm
    assert speclist[0].meta["source_id"] == "0"
    assert speclist[1].meta["source_id"] == "1"
    assert speclist[2].meta["source_id"] == "2"


data_input = ["file", "asdf", "blob"]

@pytest.fixture(params=data_input)
def data(request, roman_single):
    if request.param == "file":
        tmp = roman_single
    elif request.param == "asdf":
        tmp = asdf.open(roman_single)
    elif request.param == "blob":
        tmp = open(roman_single, 'rb')
    else:
        tmp = roman_single

    yield tmp

    try:
        tmp.close()
    except AttributeError:
        pass

def test_read_fileobj_or_asdftree(data):
    """Test Spectrum read with different data inputs"""
    s = Spectrum.read(data)
    assert isinstance(s, Spectrum)
