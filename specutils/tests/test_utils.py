from astropy import units as u
from astropy import modeling
from specutils.utils import UnitedModel


def test_united_model():
    c = modeling.models.Chebyshev1D(3)
    uc = UnitedModel(c, u.AA, u.km)

    assert uc(10*u.nm).to(u.m) == 0*u.m
