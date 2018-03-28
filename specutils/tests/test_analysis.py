import astropy.units as u
import astropy.wcs as fitswcs
import gwcs
import numpy as np
import pytest

from ..spectra.spectrum1d import Spectrum1D
from ..analysis import equivalent_width


def test_equivalent_width():
    spec = Spectrum1D(spectral_axis=np.arange(50),
                      flux=np.array([
        0.09731049,  0.04101204,  0.09726903,  0.72524453, -0.28105335,
       -0.93772961, -1.9828759 ,  0.38752423,  0.86006845, -0.08198352,
       -0.08303639,  0.18421212, -0.50724803, -0.09625829, -1.9318252 ,
        1.07973092, -0.3189966 , -1.52045995,  1.95926732,  1.71674612,
        0.28450979,  0.37737352, -1.16525665,  0.29277855, -0.37458935,
       -1.31719473, -0.31894975, -0.51095169, -0.45959643,  0.77837891,
        0.91153499,  0.13612405,  0.63433898, -0.91986964, -0.4546604 ,
       -1.09797558, -1.83641516, -0.94179757, -1.33989398, -0.66452815,
       -0.71835507, -1.39939311,  0.50070437, -1.03926682,  0.58481419,
        0.19552929, -0.7862626 ,  0.51592792, -0.95650517, -1.26917689]))

    ew = equivalent_width(spec)

    assert isinstance(ew, u.Quantity)
    assert np.allclose(ew.value, 6.8278704893358)

    ew = equivalent_width(spec[10:20])

    assert np.allclose(ew.value, 15.37809622)

