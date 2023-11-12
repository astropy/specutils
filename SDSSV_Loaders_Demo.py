r"""°°°
# New SDSS-V loaders for specutils
### Riley Thai
°°°"""
#|%%--%%| <EudVIMyXP8|CSGfqztTds>

%load_ext autoreload
%autoreload 2
#%matplotlib inline
import matplotlib.pyplot as plt
import astropy.visualization as viz
import astropy.units as u
import numpy as np

import sdss_v

# my personal fancy renderer/style files
import scienceplots
import matplotlib as mpl

#|%%--%%| <CSGfqztTds|r2SyGKMq6X>
r"""°°°
## APOGEE files
°°°"""
#|%%--%%| <r2SyGKMq6X|HUnn0AELlP>
r"""°°°
### ApStar -> Spectrum1D

Very simply loads the coadd.
°°°"""
#|%%--%%| <HUnn0AELlP|ZanVRNtl3m>

spectrum = sdss_v.load_sdss_apStar_1D(
    "/home/riley/uni/rproj/data/apStar-1.2-apo25m-2M05560393-0133549.fits")

with viz.quantity_support():
    flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
    plt.plot(spectrum.spectral_axis,np.transpose(flux))
plt.show()


print(spectrum)
#|%%--%%| <ZanVRNtl3m|7EKchx373x>
r"""°°°
## ApStar -> SpectrumList

This test file doesn't seem to have visit spectra attached with it, so I haven't verified the functionality of the List output fully.
°°°"""
#|%%--%%| <7EKchx373x|nyghrvzu7c>

spectra = sdss_v.load_sdss_apStar_list(
    "/home/riley/uni/rproj/data/apStar-1.2-apo25m-2M05560393-0133549.fits")

spectrum = spectra[0]
with viz.quantity_support():
    flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
    plt.plot(spectrum.spectral_axis,np.transpose(flux))
plt.show()

print(spectra)

#|%%--%%| <nyghrvzu7c|UyWl2KZ4CB>
r"""°°°
### apVisit -> SpectrumList

Each chip is its own Spectrum1D object. Note that no chip order information seems to be saved, so I can't get like a "chip ID" to metadata.
°°°"""
#|%%--%%| <UyWl2KZ4CB|lAPMaPnQhC>

spectra= sdss_v.load_sdss_apVisit_multi(
    "/home/riley/uni/rproj/data/apVisit-1.2-apo25m-3786-59637-275.fits")

with viz.quantity_support():
    for i,spectrum in enumerate(spectra):
        flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
        plt.plot(spectrum.spectral_axis,np.transpose(flux), label=i)
plt.legend(loc='best')
plt.show()

print(spectra)
#|%%--%%| <lAPMaPnQhC|TOEIEGGL9M>
r"""°°°
### apVisit -> Spectrum1D

Compiles all chips into a single Spectrum1D object by concatenating the 3 chips' spectra. Not sure if this is what we're after?
°°°"""
#|%%--%%| <TOEIEGGL9M|BFzUNZnaiJ>

spectrum= sdss_v.load_sdss_apVisit_1D(
    "/home/riley/uni/rproj/data/apVisit-1.2-apo25m-3786-59637-275.fits")

with viz.quantity_support():
    flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
    plt.plot(spectrum.spectral_axis,np.transpose(flux),label="{}--{}".format(spectrum.meta['mjd'],spectrum.meta['date-obs']))
plt.legend(loc='best')
plt.show()


print(spectrum)
#|%%--%%| <BFzUNZnaiJ|KBnAScSlkw>
r"""°°°
## BOSS files
°°°"""
#|%%--%%| <KBnAScSlkw|oRfGYVh8ED>
r"""°°°
### specFull -> Spectrum1D

This only loads coadd at HDU1.
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
### specFull -> SpectrumList 

Loads coadd + all exposures as Spectrum1D objects.
°°°"""
#|%%--%%| <IZ16ccvyCG|ytyxfeis77>

spectra = sdss_v.load_sdss_specFull_list("/home/riley/uni/rproj/data/spec-015252-59278-4593082715.fits")
print(type(spectra), ":", spectra)

for i,spectrum in enumerate(spectra):
    with viz.quantity_support():
        flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
        plt.plot(spectrum.spectral_axis,np.transpose(flux), label=spectrum.meta['name'],alpha=0.5,zorder=3-i)
    plt.legend(loc='best')
    plt.yscale('log')
plt.show()

#|%%--%%| <ytyxfeis77|HDZ0IOUeN7>
r"""°°°
### specLite -> Spectrum1D 

This only loads coadd at HDU1.
°°°"""
#|%%--%%| <HDZ0IOUeN7|s6hoXAakRi>

spectrum = sdss_v.load_sdss_specLite_1D("/home/riley/uni/rproj/data/spec-015252-59278-4593082715.fits")
print(type(spectrum), ":", spectrum)
with viz.quantity_support():
    flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
    plt.plot(spectrum.spectral_axis,np.transpose(flux))
    plt.yscale('log')
plt.show()


#|%%--%%| <s6hoXAakRi|as2OmK8Ksf>
r"""°°°
### specLite -> SpectrumList 

Loads coadd + all exposures as Spectrum1D objects. The DSI says it can have exposures, but most likely it won't/shouldn't, so this is a low priority case.
°°°"""
#|%%--%%| <as2OmK8Ksf|KySoNu77pr>

spectra = sdss_v.load_sdss_specLite_list("/home/riley/uni/rproj/data/spec-015252-59278-4593082715.fits")
print(type(spectra), ":", spectra)

for i,spectrum in enumerate(spectra):
    with viz.quantity_support():
        flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
        plt.plot(spectrum.spectral_axis,np.transpose(flux), label=spectrum.meta['name'],alpha=0.5,zorder=3-i)
    plt.legend(loc='best')
    plt.yscale('log')
plt.show()

#|%%--%%| <KySoNu77pr|MyuQK6fz8Q>
r"""°°°
### Custom Coadds (allEpoch, etc) -> specFull cases

I like this BOSS team's organization. Everything is formatted the same (very happy).
°°°"""
#|%%--%%| <MyuQK6fz8Q|MlmlHgCOYQ>

spectra = sdss_v.load_sdss_specFull_list("/home/riley/uni/rproj/data/spec-allepoch-60130-27021597775058535.fits")
print(type(spectra), ":", spectra)

for i,spectrum in enumerate(spectra):
    with viz.quantity_support():
        flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
        plt.plot(spectrum.spectral_axis,np.transpose(flux), label=spectrum.meta['name'],alpha=0.5,zorder=3-i)
    plt.legend(loc='best')
    plt.yscale('log')
plt.show()

#|%%--%%| <MlmlHgCOYQ|IFjhstzMoW>

spectra = sdss_v.load_sdss_specFull_list("/home/riley/uni/rproj/data/spec-112359-60118-7612895412.fits")
print(type(spectra), ":", spectra)

for i,spectrum in enumerate(spectra):
    with viz.quantity_support():
        flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
        plt.plot(spectrum.spectral_axis,np.transpose(flux), label=spectrum.meta['name'],alpha=0.5,zorder=3-i)
    plt.legend(loc='best')
    plt.yscale('log')
plt.show()

#|%%--%%| <IFjhstzMoW|MdP1J8W9oG>
r"""°°°
## ASTRA files (MWM)
°°°"""
#|%%--%%| <MdP1J8W9oG|XX2R7Sy1CQ>
r"""°°°
### mwmVisit -> Spectrum1D

Spectrum1D of resampled rest-frame data from a SINGLE type of spectrograph (based on the HDU specified).

Within the meta, there are arrays of size n, where n is the number of visits.
°°°"""
#|%%--%%| <XX2R7Sy1CQ|XptfhNXZsC>

spectrum = sdss_v.load_sdss_mwmVisit_1d("/home/riley/uni/rproj/data/mwmVisit-0.5.0-70350000.fits",3)
print(type(spectrum), ":", spectrum)
spectral_axis = spectrum.spectral_axis
flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
for i in range(len(flux)):
    with viz.quantity_support():
        plt.plot(spectral_axis,np.transpose(flux[i]),label=spectrum.meta["date"][i])
plt.title(spectrum.meta['name'])
plt.legend(loc='best')
plt.show()

#|%%--%%| <XptfhNXZsC|bIMJg22G3T>
r"""°°°
### mwmVisit -> SpectrumList[Spectrum1D]

A list of spectra, where each Spectrum1D object contains all the flux and meta for each visit.
°°°"""
#|%%--%%| <bIMJg22G3T|zt0s14z7ro>

spectra = sdss_v.load_sdss_mwmVisit_list("/home/riley/uni/rproj/data/mwmVisit-0.5.0-70350000.fits")

fig, axes = plt.subplots(nrows=2,ncols=2,layout='constrained')
axes = axes.flatten()
locations = ['BOSS/APO','BOSS/LCO',"APOGEE/APO","APOGEE/LCO"]
for q, spectrum in enumerate(spectra):
    if spectrum is None:
        axes[q].set_title(locations[q])
        axes[q].text(0.5,0.5,'No Spectra',horizontalalignment='center',verticalalignment='center')
        continue
    spectral_axis = spectrum.spectral_axis
    flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
    for i in range(len(flux)):
        with viz.quantity_support():
            axes[q].plot(spectral_axis,np.transpose(flux[i]),label=spectrum.meta["date"][i])
    axes[q].legend(loc='best')
    axes[q].set_title(spectra[q].meta['type'])
fig.show()
#|%%--%%| <zt0s14z7ro|rEGLOBO3Q9>
r"""°°°
### mwmStar -> Spectrum1D

Spectrum1D, only loading 1 of the specified HDU's
°°°"""
#|%%--%%| <rEGLOBO3Q9|fdJCv1YhzM>

spectrum = sdss_v.load_sdss_mwmStar_1d("/home/riley/uni/rproj/data/mwmStar-0.5.0-103020000.fits",4)
print(type(spectrum), ":", spectrum)
spectral_axis = spectrum.spectral_axis
flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
with viz.quantity_support():
    plt.plot(spectral_axis,np.transpose(flux),label=spectrum.meta["mjd"],color='black')
plt.legend(loc='best')
plt.title(spectrum.meta['data'])

plt.show()


#|%%--%%| <fdJCv1YhzM|Q64T9opCSO>
r"""°°°
### mwmStar -> SpectrumList[Spectrum1D]

SpectrumList of Spectrum1D, like mwmVisit. Lists all valid HDU's.
°°°"""
#|%%--%%| <Q64T9opCSO|yMudsqYfdk>

spectra = sdss_v.load_sdss_mwmStar_list("/home/riley/uni/rproj/data/mwmStar-0.5.0-103020000.fits")

fig, axes = plt.subplots(nrows=2,ncols=2,layout='constrained')
axes = axes.flatten()
locations = ['BOSS/APO','BOSS/LCO',"APOGEE/APO","APOGEE/LCO"]
for q, spectrum in enumerate(spectra):
    if spectrum is None:
        axes[q].set_title(locations[q])
        axes[q].text(0.5,0.5,'No Spectra',horizontalalignment='center',verticalalignment='center')
        continue
    spectral_axis = spectrum.spectral_axis
    flux = spectrum.flux.to("erg / (Angstrom cm2 s)")
    with viz.quantity_support():
        axes[q].plot(spectral_axis,np.transpose(flux),label=spectrum.meta["mjd"],color='black')
        axes[q].legend(loc='best')
    axes[q].set_title(spectra[q].meta['type'])
fig.show()

#|%%--%%| <yMudsqYfdk|mpdYgzpxdn>
r"""°°°
I don't have a specLite file to test on, but I'd presume I can just use the specFull 1D case as the specLite loader.

In any case, this is mainly to see if I'm loading these files as we'd expect, and just to update. I'm going to start work on the Solara stuff now.

**Current other todos:**

- Add documentation to each loader function.

- Ensure code meets standards.

- Ask someone for fake data for unit tests.

- Ensure WCS information is provided during loader in standard FITS format (read in docs something about SDSS using weird WCS format).

- Work on the identifier function for the main data loader within specutils.

    - I'll begin working on this identifier once we're sure that this is all correct. 
    

If you could tell me about the file naming conventions used for these data products, that would be great to help me get started on the identifier. 

Also, I still need the other data types (custom coadds, specLite, etc) for verification and testing.

**Some of my notes and questions:**

- How do we distinctly identify visits to a user within its metadata?

  - Some of these files contain say 5 APOGEE spectra -- how do we wish to identify it?

    - Currently, I've left these as just `DATE-OBS` date for mwm things, and occassionaly `MJD5` if I found that the column existed, but this gets messy once we start considering coadds.

  - What should the identifiers be for BOSS files, or are we leaving it as the MJD + OBJECT + FIBER combination as Joel said?

  - I presume in future, `ASTRA` files will have some sort of output identifier aside from CHPC unid. Any pre-emptive ideas on what this will be coming up?

- What is the flux unit (BUNIT) info stored as in MWM Files? Where is it stored (if it is)?
°°°"""
#|%%--%%| <mpdYgzpxdn|lkl6AQvnpl>
r"""°°°
Other testing -- ignore
°°°"""
#|%%--%%| <lkl6AQvnpl|ybE9h22peP>
from astropy.io import fits
image = fits.open("/home/riley/uni/rproj/data/mwmStar-0.3.0-27021597838600303.fits")

