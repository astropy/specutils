# Specutils

Specutils is an Astropy affiliated package with the goal of providing a shared 
set of Python representations of astronomical spectra and basic tools to 
operate on these spectra. The effort is also meant to be a "hub", helping to 
unite the Python astronomical spectroscopy community around shared effort, 
much as Astropy is meant to for the wider astronomy Python ecosystem.

Note that Specutils is not intended as an all-in-one spectroscopic analysis or 
reduction tool. While it provides some basic analysis (following the Python 
philosophy of "batteries included"), it is also meant to facilitate connecting 
together disparate reduction pipelines and analysis tools through shared data 
representations.


## Installation

There are several ways to install the package, the most direct is using `pip`

```bash
$ pip install git+https://github.com/nmearl/specutils
```

Otherwise, you may simply clone the repo and install the package

```bash
$ git clone https://github.com/nmearl/specutils
$ cd specutils
$ python setup.py install
```

or, if you'd like to easily be able to `pip uninstall`, use the following 
commands

```bash
$ git clone https://github.com/nmearl/specutils
$ cd specutils
$ pip install .
```

## Quickstart

Defining a spectrum is straightforward

```python
from astropy.units import Quantity
from specutils.spectra import Spectrum1D
import numpy as np

# Using Astropy `Quantity`s
spec = Spectrum1D(flux=Quantity(np.random.sample(100), "erg/Angstrom/cm2/s"),
                  dispersion=Quantity(np.arange(100), "Angstrom"))

# Without `Quantity`s
spec = Spectrum1D(flux=np.random.sample(100), unit="erg/Angstrom/cm2/s",
                  dispersion=np.arange(100), disp_unit="Angstrom")
```

Converting to different units

```python
spec.to_flux("Jy")

# From Angstrom to Hz
spec.to_dispersion("Hz", rest=Quantity(202, "Angstrom"))

# From Hz to cm/s
spec.to_dispersion("cm/s", rest=Quantity(1.48412108e+16, "Hz"))

# From cm/s to Angstrom
spec.to_dispersion("Angstrom", rest=Quantity(1.48412108e+16, "Hz"))
```

Basic analysis usage

```python
from specutils.analysis import equivalent_width

print(equivalent_width(spec))
```

### Custom loaders

Define a custom loader in a separate python file and place the file in your
`~/.specutils` directory. Upon importing `specutils`, the loader will be added
to the registry.

```python
# ~/.specutils/my_custom_loader.py
import os
import six

from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.table import Table
from astropy.units import Unit
from astropy.wcs import WCS

from specutils.io.registers import data_loader
from specutils.spectra import Spectrum1D

# Define an optional identifier. If made specific enough, this circumvents the
# need to add `format="my-format"` in the `Spectrum1D.read` call.
def identify_generic_fits(origin, *args, **kwargs):
    return (isinstance(args[0], six.string_types) and
            os.path.splitext(args[0].lower())[1] == '.fits')


@data_loader("my-format", identifier=identify_generic_fits)
def generic_fits(file_name, **kwargs):
    name = os.path.basename(file_name.rstrip(os.sep)).rsplit('.', 1)[0]

    with fits.open(file_name, **kwargs) as hdulist:
        header = hdulist[0].header

        tab = Table.read(file_name)

        meta = {'header': header}
        wcs = WCS(hdulist[0].header)
        uncertainty = StdDevUncertainty(tab["err"])
        data = tab["flux"] * Unit("Jy")

    return Spectrum1D(flux=data, wcs=wcs, uncertainty=uncertainty, meta=meta)
```

Using your custom loader:

```python
from specutils import Spectrum1D

spec = Spectrum1D("path/to/data", format="my-format")
```