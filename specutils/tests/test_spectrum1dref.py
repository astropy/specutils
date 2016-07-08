from specutils.core.generic import Spectrum1DRef
import numpy as np
import astropy.units as u
from astropy.wcs import WCS


def test_spectrum1dref_from_array():
    """
    Test the ability to create a spectrum object from data and dispersion
    arrays with optional units.
    """
    data = np.random.sample(100)
    disp = np.arange(100)
    spectrum = Spectrum1DRef.from_array(data, dispersion=disp,
                                        unit=u.Unit("Jy"),
                                        dispersion_unit=u.Unit("Angstrom"))

    # Compare units
    assert spectrum.unit == u.Unit("Jy")
    assert spectrum.dispersion_unit == u.Unit("Angstrom")

    # Test data
    assert (spectrum.data == data).all()
    assert (spectrum.dispersion == disp).all()


def test_spectrum1dref_wcs():
    """
    Test the ability to extract information from wcs.
    """
    wcs = WCS({"CUNIT1": "Angstrom", "CRVAL1": 200.0, "CDELT1": 1.0})
    data = np.random.sample(100)
    spectrum = Spectrum1DRef(data=data, wcs=wcs,
                             unit=u.Unit("erg/cm^2/s/Angstrom"))

    # Check the dispersion constructed from the WCS
    assert spectrum.dispersion[0] == 201.0

    # Check the unit information taken from the WCS
    assert spectrum.unit == u.Unit("erg/cm^2/s/Angstrom")
    assert spectrum.dispersion_unit == u.Unit("Angstrom")


def test_spectrum1dref_arithmetic():
    """
    Test arithmetic functionality.
    """
    data1 = np.random.sample(100)
    data2 = np.random.sample(100)

    spectrum1 = Spectrum1DRef(data1, unit=u.Unit("Jy/cm/s"))
    spectrum2 = Spectrum1DRef(data2, unit=u.Unit("erg/cm^2/s/Angstrom"))

    assert (spectrum1 + spectrum2).unit == spectrum1.unit
    assert (spectrum1 + spectrum2).data[0] == \
           (data1 * u.Unit("Jy/cm/s") + data2 * u.Unit("erg/cm^2/s/Angstrom")).value[0]


def test_spectrum1dref_slicing():
    """
    Test ability to slice the Spectrum1DRef data object.
    """
    data = np.random.sample(100)
    spectrum = Spectrum1DRef(data)

    assert spectrum[0].data == data[0]


def test_spectrum1dref_copy():
    """
    Test the ability to duplicate a Spectrum1DRef object.
    """
    wcs = WCS({"CUNIT1": "Angstrom", "CRVAL1": 200.0, "CDELT1": 1.0})
    data = np.random.sample(100)
    disp = np.arange(100)
    spectrum = Spectrum1DRef.from_array(data, dispersion=disp,
                                        unit=u.Unit("Jy"),
                                        dispersion_unit=u.Unit("Angstrom"),
                                        wcs=wcs, meta={'author': 'me'})

    spectrum_copy = Spectrum1DRef.copy(spectrum)

    # Check data
    assert (spectrum.data == spectrum_copy.data).all()
    assert (spectrum.dispersion == spectrum_copy.dispersion).all()

    # Check unit, wcs, meta
    assert spectrum.unit == spectrum_copy.unit
    assert spectrum.wcs.wcs.cunit[0] == spectrum_copy.wcs.wcs.cunit[0]
    assert spectrum.wcs.wcs.crval[0] == spectrum_copy.wcs.wcs.crval[0]
    assert spectrum.wcs.wcs.cdelt[0] == spectrum_copy.wcs.wcs.cdelt[0]
    assert spectrum.meta == spectrum_copy.meta