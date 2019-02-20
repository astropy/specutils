# -*- coding: utf-8 -*-

import pytest
import numpy as np
from astropy import constants as const

from ..manipulation.shift import gravitational_redshift

def test_gravitational_redshift():
    assert np.isclose(gravitational_redshift(solar=True).value, 2.1119411024450585e-06)
    zsolar = gravitational_redshift(solar=False, obj_mass=const.M_sun,
                                    obj_radius=const.R_sun, obj_distance=const.au)
    assert gravitational_redshift(solar=True) == zsolar
    with pytest.raises(AssertionError):
        gravitational_redshift(solar=False, obj_mass=1)
