# tests for FITS WCS IO

from astropy.io import fits

simple_linear_wcs_header = fits.Header([('ctype1', 'LINEAR'), ('crval1', 6000.), ('crpix1', 1), ('cd1_1', 1.0),
                                        ('cdelt1', 1.0)])
