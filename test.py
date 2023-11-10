r"""°°°
Hello!
°°°"""
#|%%--%%| <EudVIMyXP8|CSGfqztTds>

%load_ext autoreload
%autoreload 2
%matplotlib inline
import matplotlib.pyplot as plt
import scienceplots
import astropy.visualization as viz
import astropy.units as u
import numpy as np

import sdss_v

import matplotlib as mpl
mpl.use("module://matplotlib-backend-kitty")

#|%%--%%| <CSGfqztTds|HUnn0AELlP>
r"""°°°
ApStar.fits -> Spectrum1D
°°°"""
#|%%--%%| <HUnn0AELlP|1VXcbmCQDE>

spectrum = sdss_v.load_sdss_apStar_1D(
    "/home/riley/uni/rproj/data/apStar-1.2-apo25m-2M05560393-0133549.fits")
print(spectrum)
print(type(spectrum))
#|%%--%%| <1VXcbmCQDE|FFfDu7vYi7>

with viz.quantity_support():
    flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
    plt.plot(spectrum.spectral_axis,np.transpose(flux))
plt.show()

#|%%--%%| <FFfDu7vYi7|7EKchx373x>
r"""°°°
ApStar -> SpectrumList 

For some reason, this test file doesn't seem to have visit spectra.
°°°"""
#|%%--%%| <7EKchx373x|Fru1CCLRFU>

spectra = sdss_v.load_sdss_apStar_list(
    "/home/riley/uni/rproj/data/apStar-1.2-apo25m-2M05560393-0133549.fits")
print(spectra)
print(type(spectra))
print(type(spectra[0]))
#|%%--%%| <Fru1CCLRFU|POT0OZmbyH>

spectrum = spectra[0]
with viz.quantity_support():
    flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
    plt.plot(spectrum.spectral_axis,np.transpose(flux))
plt.show()

#|%%--%%| <POT0OZmbyH|UyWl2KZ4CB>
r"""°°°
apVisit -> SpectrumList
°°°"""
#|%%--%%| <UyWl2KZ4CB|R4cE41Xkjd>

spectra= sdss_v.load_sdss_apVisit_multi(
    "/home/riley/uni/rproj/data/apVisit-1.2-apo25m-3786-59637-275.fits")
print(spectra)
print(type(spectra))

#|%%--%%| <R4cE41Xkjd|45zqQFCHg3>

with viz.quantity_support():
    for spectrum in spectra:
        flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
        plt.plot(spectrum.spectral_axis,np.transpose(flux))
plt.show()
#|%%--%%| <45zqQFCHg3|TOEIEGGL9M>
r"""°°°
apVisit -> Spectrum1D

compiles all chips into a single Spectrum1D object by concatenating the 3 chips' spectra. Not sure if this is what we're after?
°°°"""
#|%%--%%| <TOEIEGGL9M|ApSShh1BtL>

spectrum= sdss_v.load_sdss_apVisit_1D(
    "/home/riley/uni/rproj/data/apVisit-1.2-apo25m-3786-59637-275.fits")
print(spectrum)
print(type(spectrum))

#|%%--%%| <ApSShh1BtL|PkHUQw0tM3>

with viz.quantity_support():
    flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
    plt.plot(spectrum.spectral_axis,np.transpose(flux))
plt.show()

#|%%--%%| <PkHUQw0tM3|oRfGYVh8ED>
r"""°°°
specFull -> Spectrum1D (only loads coadd at HDU1)
°°°"""
#|%%--%%| <oRfGYVh8ED|l20TC18Pgi>

spectrum = sdss_v.load_sdss_specFull_1D("/home/riley/uni/rproj/data/spec-015252-59278-4593082715.fits")
print(type(spectrum), ":", spectrum)
with viz.quantity_support():
    flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
    plt.plot(spectrum.spectral_axis,np.transpose(flux))
    plt.yscale('log')
plt.show()

#|%%--%%| <l20TC18Pgi|IZ16ccvyCG>
r"""°°°
specFull -> SpectrumList (coadd + all exposures)
°°°"""
#|%%--%%| <IZ16ccvyCG|g9HndzCSN5>

spectra = sdss_v.load_sdss_specFull_list("/home/riley/uni/rproj/data/spec-015252-59278-4593082715.fits")
print(type(spectra), ":", spectra)

#|%%--%%| <g9HndzCSN5|u8gNxVShP0>

for i,spectrum in enumerate(spectra):
    with viz.quantity_support():
        flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
        plt.plot(spectrum.spectral_axis,np.transpose(flux), label=spectrum.meta['name'],alpha=0.5,zorder=3-i)
    plt.legend(loc='best')
    plt.yscale('log')
plt.show()

#|%%--%%| <u8gNxVShP0|lkl6AQvnpl>
r"""°°°
Other testing -- ignore
°°°"""
#|%%--%%| <lkl6AQvnpl|ybE9h22peP>
from astropy.io import fits
image = fits.open("/home/riley/uni/rproj/data/spec-015252-59278-4593082715.fits")
image2 = fits.open("/home/riley/uni/rproj/data/apStar-1.2-apo25m-2M05560393-0133549.fits")

