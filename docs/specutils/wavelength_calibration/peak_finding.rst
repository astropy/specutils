Peak Finding
------------

.. warning::
    As this is my first attempt at using RST, this document will very much be
    a work in progress. I'll figure things out as I go along.

=================
Import statements
=================

>>> import matplotlib.pyplot as plt
>>> import pylab

================
Loading the data
================

There is more than one way to load FITS data into Python. Here is one method
that involves ds9; an advantage of this method is that the data can easily be
visualized as it is being manipulated.

First of all, instructions for downloading and installing the pyds9 module
can be found here: http://hea-www.harvard.edu/RD/ds9/pyds9/

Start Python, then execute the following commands:

>>> from ds9 import *
>>> d = ds9()

To read FITS data or a raw array from ds9 into pyfits, use the ‘get_pyfits’ method. It takes no args and returns an hdu list:

>>> d.set("file ./Downloads/arc.fits")
>>> hdul = d.get_pyfits()
>>> data = hdul[0].data

Collapse the data:

>>> arcsum = data.sum(axis=0)

Find the peaks:

>>> import numpy as np
>>> from scipy import signal
>>> peakind = signal.find_peaks_cwt(arcsum, np.arange(1,10))

=================
Plot the spectrum
=================

>>> import matplotlib.pyplot as plt
>>> import pylab
>>> plt.plot(arcsum)
>>> for i in peakind:
>>>     plt.vlines(i, 0, arcsum[i], color='red')
>>> pylab.show()

=================
With another file
=================

Let's put it all together with another example.

>>> import matplotlib.pyplot as plt
>>> import pylab
>>> from ds9 import *
>>> import numpy as np
>>> from scipy import signal
>>> import matplotlib.pyplot as plt
>>> import pylab

>>> d = ds9()
>>> d.set("file ./Downloads/ftlrs360254.fits")
>>> hdul2 = d.get_pyfits()

>>> data2 = hdul2[0].data
>>> arcsum2 = data2.sum(axis=0)

>>> peakind2 = signal.find_peaks_cwt(arcsum2, np.arange(1,10))

>>> plt.plot(arcsum2)
>>> for i in peakind2:
>>>     plt.vlines(i,0,arcsum2[i],color='red')
>>> pylab.show()
