r"""°°°
Hello!
°°°"""
#|%%--%%| <EudVIMyXP8|T3cto1hBOd>

%load_ext autoreload
%autoreload 2
%matplotlib inline
import matplotlib.pyplot as plt
import scienceplots
import astropy.visualization as viz
from astropy.io import fits
import numpy as np

import sdss_v

import matplotlib as mpl
mpl.use("module://matplotlib-backend-kitty")

#|%%--%%| <T3cto1hBOd|HUnn0AELlP>
r"""°°°
ApStar.fits -> Spectrum1D
°°°"""
# |%%--%%| <HUnn0AELlP|HOCivwpkH2>

spectrum = sdss_v.load_sdss_apStar(
    "/home/riley/Downloads/apStar-1.2-apo25m-2M05560393-0133549.fits")
print(spectrum)
print(type(spectrum))
#|%%--%%| <HOCivwpkH2|FFfDu7vYi7>

with viz.quantity_support():
    plt.plot(spectrum.spectral_axis,np.transpose(spectrum.flux))
plt.show()

#|%%--%%| <FFfDu7vYi7|7EKchx373x>
r"""°°°
ApStar -> SpectrumList 

For some reason, this test file doesn't seem to have visit spectra.
°°°"""
#|%%--%%| <7EKchx373x|Fru1CCLRFU>

spectra = sdss_v.load_sdss_apStar_list(
    "/home/riley/Downloads/apStar-1.2-apo25m-2M05560393-0133549.fits")
print(spectra)
print(type(spectra))
print(type(spectra[0]))
#|%%--%%| <Fru1CCLRFU|POT0OZmbyH>

spectrum = spectra[0]
with viz.quantity_support():
    plt.plot(spectrum.spectral_axis,np.transpose(spectrum.flux))
plt.show()

#|%%--%%| <POT0OZmbyH|UyWl2KZ4CB>
r"""°°°
apVisit -> SpectrumList
°°°"""
#|%%--%%| <UyWl2KZ4CB|R4cE41Xkjd>

spectra= sdss_v.load_sdss_apVisit_multi(
    "/home/riley/Downloads/apVisit-1.2-apo25m-3786-59637-275.fits")
print(spectra)
print(type(spectra))

#|%%--%%| <R4cE41Xkjd|45zqQFCHg3>

with viz.quantity_support():
    for spectrum in spectra:
        plt.plot(spectrum.spectral_axis,np.transpose(spectrum.flux))
plt.show()
#|%%--%%| <45zqQFCHg3|TOEIEGGL9M>
r"""°°°
apVisit -> Spectrum1D

compiles all chips into a single spectra by concatenating the 3 chips. Not sure if this is what we're after?
°°°"""
#|%%--%%| <TOEIEGGL9M|ApSShh1BtL>

spectrum= sdss_v.load_sdss_apVisit(
    "/home/riley/Downloads/apVisit-1.2-apo25m-3786-59637-275.fits")
print(spectrum)
print(type(spectrum))

#|%%--%%| <ApSShh1BtL|PkHUQw0tM3>

with viz.quantity_support():
    plt.plot(spectrum.spectral_axis,np.transpose(spectrum.flux))
plt.show()

#|%%--%%| <PkHUQw0tM3|lkl6AQvnpl>
r"""°°°
Other testing -- ignore
°°°"""
#|%%--%%| <lkl6AQvnpl|ybE9h22peP>

image = fits.open("/home/riley/Downloads/apVisit-1.2-apo25m-3786-59637-275.fits")

