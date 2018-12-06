from astropy import units as u
from astropy import modeling
from specutils.utils import QuantityModel


def test_quantity_model():
    c = modeling.models.Chebyshev1D(3)
    uc = QuantityModel(c, u.AA, u.km)

    assert uc(10*u.nm).to(u.m) == 0*u.m
