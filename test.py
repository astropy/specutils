r"""°°°
Hello!
°°°"""
#|%%--%%| <EudVIMyXP8|4kasgaDERv>

%load_ext autoreload
%autoreload 2
%matplotlib inline
import matplotlib.pyplot as plt
import matplotlib as mpl
import scienceplots
import astropy.visualization as viz
import numpy as np

#|%%--%%| <4kasgaDERv|HUnn0AELlP>
r"""°°°
ApStar.fits -> Spectrum1D
°°°"""
# |%%--%%| <HUnn0AELlP|HOCivwpkH2>
from sdss_v import data_loader, load_sdss_apStar

spectra = load_sdss_apStar(
    "/home/riley/Downloads/apStar-1.2-apo25m-2M05560393-0133549.fits")
print(spectra)
print(type(spectra))
#|%%--%%| <HOCivwpkH2|FFfDu7vYi7>

mpl.use('module://matplotlib-backend-kitty')
plt.plot(spectra.spectral_axis,np.transpose(spectra.flux))
plt.show()

#|%%--%%| <FFfDu7vYi7|7EKchx373x>
r"""°°°
ApStar -> SpectrumCollection

For some reason, this file doesn't seem to have visit spectra.
°°°"""
#|%%--%%| <7EKchx373x|Fru1CCLRFU>
from sdss_v import load_sdss_apStar_list

spectra = load_sdss_apStar_list(
    "/home/riley/Downloads/apStar-1.2-apo25m-2M05560393-0133549.fits")
print(spectra)
print(type(spectra))
print(type(spectra[0]))
#|%%--%%| <Fru1CCLRFU|POT0OZmbyH>

spectrum = spectra[0]
plt.plot(spectrum.spectral_axis,np.transpose(spectrum.flux))
plt.show()

#|%%--%%| <POT0OZmbyH|UyWl2KZ4CB>
r"""°°°
apVisit -> Spectrum1D or SpectrumList
°°°"""
#|%%--%%| <UyWl2KZ4CB|R4cE41Xkjd>

from sdss_v import load_sdss_apVisit

spectrum= load_sdss_apVisit(
    "/home/riley/Downloads/apVisit-1.2-apo25m-3786-59637-275.fits")
print(spectrum)
print(type(spectrum))

#|%%--%%| <R4cE41Xkjd|P9bX4Azf27>

plt.plot(spectrum.spectral_axis,np.transpose(spectrum.flux))
plt.show()

