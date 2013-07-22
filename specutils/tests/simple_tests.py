from specutils import Spectrum1D
from astropy import units
import numpy as np
dispersion = np.arange(4000, 5000, 0.12)
flux = np.random.randn(len(dispersion))
mySpectrum = Spectrum1D.from_array(dispersion,
                                   flux,
                                   dispersion_unit=units.m)

hBeta = mySpectrum.slice_dispersion(4851.0, 4871.0)
hBeta
