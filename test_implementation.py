r"""°°°
# Specutils implementation test notebook
### Riley Thai
°°°"""
# |%%--%%| <9aEgDHplIU|q1E77jiDZi>
from specutils import Spectrum1D, SpectrumList

from jdaviz import Specviz

Spectrum1D.read.list_formats()
SpectrumList.read.list_formats()

dir = "/home/riley/uni/rproj/data/"

# |%%--%%| <q1E77jiDZi|qUw0OFqWcl>
r"""°°°
## Spectrum1D loaders
°°°"""
# |%%--%%| <qUw0OFqWcl|vbpxGYnPFF>

specviz = Specviz()
names = [
    "mwmVisit", "mwmStar", "apStar", "apVisit", "spec1", "spec2",
    "custom coadd"
]
specs_1D = list()
specs_1D.append(Spectrum1D.read(f"{dir}mwmVisit-0.5.0-70350000.fits", hdu=3))
specs_1D.append(Spectrum1D.read(f"{dir}mwmStar-0.5.0-103020000.fits", hdu=4))
specs_1D.append(
    Spectrum1D.read(f"{dir}apStar-1.2-apo25m-2M05560393-0133549.fits"))
specs_1D.append(
    Spectrum1D.read(f"{dir}apVisit-1.2-apo25m-3786-59637-275.fits"))
specs_1D.append(
    Spectrum1D.read(f"{dir}spec-015252-59278-4593082715.fits", hdu=1))
specs_1D.append(
    Spectrum1D.read(f"{dir}spec-112359-60118-7612895412.fits", hdu=1))
specs_1D.append(
    Spectrum1D.read(f"{dir}spec-allepoch-60130-27021597775058535.fits", hdu=1))

# for i in range(7):
#    specviz.load_data(specs_1D[i], data_label=names[i])

specviz.load_data(specs_1D[3], data_label="test")

# specviz.show()

# |%%--%%| <vbpxGYnPFF|Zccee01Dam>
r"""°°°
## SpectrumList loaders
°°°"""
# |%%--%%| <Zccee01Dam|xXJL3a0cig>

specs_multi = list()
specviz_multi = Specviz()
names = ["mwmVisit", "mwmStar", "apVisit", "spec1", "spec2", "custom coadd"]
specs_multi.append(
    SpectrumList.read(f"{dir}mwmVisit-0.5.0-70350000.fits",
                      format="SDSS-V mwm multi"))
specs_multi.append(
    SpectrumList.read(f"{dir}mwmStar-0.5.0-103020000.fits",
                      format="SDSS-V mwm multi"))
specs_multi.append(
    SpectrumList.read(f"{dir}apVisit-1.2-apo25m-3786-59637-275.fits",
                      format="SDSS-V apVisit multi"))
specs_multi.append(
    SpectrumList.read(f"{dir}spec-015252-59278-4593082715.fits",
                      format="SDSS-V spec multi"))

specs_multi.append(
    SpectrumList.read(f"{dir}spec-112359-60118-7612895412.fits",
                      format="SDSS-V spec multi"))

specs_multi.append(
    SpectrumList.read(f"{dir}spec-allepoch-60130-27021597775058535.fits",
                      format="SDSS-V spec multi"))

for i in range(6):
    specviz_multi.load_data(specs_multi[i], data_label=names[i])

# specviz_multi.show()
# |%%--%%| <xXJL3a0cig|hDooH3VxTr>

specs_multi.append(
    SpectrumList.read(f"{dir}apStar-1.2-apo25m-2M05560393-0133549.fits",
                      format="SDSS-V apStar multi"))
