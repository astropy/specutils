from astropy.io.fits import Header
import pytest

from ..wcs_fits import _read_non_linear_iraf_wcs


MULTISPEC_LINEAR_HEADER = Header(
    {
        "SIMPLE": True,
        "BITPIX": -32,
        "NAXIS": 2,
        "NAXIS1": 256,
        "NAXIS2": 3,
        "WATO_001": "system-multispec",
        "WAT1_001": "wtype=multispec label-Wavelength units-Angstroms",
        "WAT2_001": 'wtype=multispec spec1 = "1 113 0 4955.44287109375 0.0568952970206737',
        "WAT2_002": '5 256 0. 23.22 31.27" spec2 = "2 112 0 4999.081054687501 0.063871018',
        "WAT2_003": '58854293 256 0. 46.09 58.44" spec3 = "3 111 0 5043.505859375 0.07096',
        "WAT2_004": '928358078002 256 0. 69.28 77.89"',
        "WCSDIM": 2,
        "DC-FLAG": 0,
        "CTYPE1": "MULTISPE",
        "LTM1_1": 1.0,
        "CD1_1": 1.0,
        "CTYPE2": "MULTISPE",
        "LTM2_2": 1.0,
        "CD2_2": 1.0,
    }
)


def test_multispec_linear_wcs():
    """Test that the multispec linear WCS is read correctly.

    To keep this test fast and short, do not use a full fits file,
    but just pass in the header.

    The example is taken from Fig 3 in the IFAF WCS spec
    https://noirlab.edu/science/sites/default/files/media/archives/documents/scidoc3040.pdf
    """
    out = _read_non_linear_iraf_wcs(MULTISPEC_LINEAR_HEADER, 2)
    assert out[0][0:3] == pytest.approx([4955.44287109, 4955.49976639, 4955.55666169])


MULTISPEC_LOGLINEAR_HEADER = Header(
    {
        "SIMPLE": True,
        "BITPIX": -32,
        "NAXIS": 2,
        "NAXIS1": 256,
        "NAXIS2": 3,
        "WATO_001": "system-multispec",
        "WAT1_001": "wtype=multispec label-Wavelength units-Angstroms",
        "WAT2_001": 'wtype=multispec spec1 = "1 113 1 4.000 0.0100 256 0. 23.22 31.27" ',
        "WAT2_002": 'spec2 = "2 112 1 4.010 0.0100 256 0. 46.09 58.44"',
        "WCSDIM": 2,
        "DC-FLAG": 0,
        "CTYPE1": "MULTISPE",
        "LTM1_1": 1.0,
        "CD1_1": 1.0,
        "CTYPE2": "MULTISPE",
        "LTM2_2": 1.0,
        "CD2_2": 1.0,
    }
)


def test_multispec_loglinear_wcs():
    """Test that the multispec loglinear WCS is read correctly.

    To keep this test fast and short, do not use a full fits file,
    but just pass in the header.
    """
    out = _read_non_linear_iraf_wcs(MULTISPEC_LOGLINEAR_HEADER, 2)
    assert out[0][0:3] == pytest.approx([10000., 10232.92992281, 10471.28548051])
    assert out[1][0:3] == pytest.approx([10232.92992281, 10471.28548051, 10715.19305238])
