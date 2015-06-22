from pylab import *
from specutils.io import read_fits
myspec = read_fits.read_fits_spectrum1d('../../specutils/io/tests/files/UVES.fits', dispersion_unit='angstrom')
plot(myspec.wavelength, myspec.flux)
show()
