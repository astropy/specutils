Peak Finding
------------

Peak finding is a task that I will be writing documentation about.

.. warning::
    As this is my first attempt at using RST, this document will very much be
    a work in progress. I'll figure things out as I go along.

=================
Import statements
=================

Don't forget these - I need numpy, matplotlib/pylab, and signals at least.

================
Loading the data
================

There is more than one way to load FITS data into Python. Here is one method
that involves ds9; an advantage of this method is that the data can easily be
visualized as it is being manipulated.

To read FITS data or a raw array from ds9 into pyfits, use the ‘get_pyfits’ method. It takes no args and returns an hdu list:

>>> d.set("file ./Downloads/arc.fits")
>>> hdul = d.get_pyfits()
>>> hdul.info()
Filename: StringIO.StringIO
No.    Name         Type      Cards   Dimensions   Format
0    PRIMARY     PrimaryHDU      24  (1024, 1024)  float32
>>> data = hdul[0].data
>>> data.shape
(1024, 1024)

Collapse the data by doing this:

>>> arcsum = data.sum(axis=0)

Find the peaks:

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

>>> d.set("file ./Downloads/ftlrs360254.fits")
>>> hdul2 = d.get_pyfits()
>>> data2 = hdul2[0].data
>>> arcsum2 = data2.sum(axis=0)
>>> peakind2 = signal.find_peaks_cwt(arcsum2, np.arange(1,10))
>>> arcsum
>>> len(arcsum)
>>> len(arcsum2)
>>> data2
>>> data
>>> plt.plot(arcsum2)
>>> for i in peakind2:
>>>     plt.vlines(i,0,arcsum2[i],color='red')
>>> pylab.show()

