###########################
Working with Spectral Cubes
###########################

Spectral cubes can be read in using the
`SpectralCube <https://spectral-cube.readthedocs.io/en/latest/>`_ package then
collapsed using various methods to yield a Spectrum1D object.  A specific
example of this workflow is demonstrated here.

Note that the workflow described here is for spectral cubes that are rectified
such that one of the axes is entirely spectral and all the spaxels have the same
``spectral_axis`` values - i.e., case 2 in
:ref:`specutils-representation-overview`. For less-rectified cubes,
pre-processing steps (not addressed by specutils at the time of this writing)
will be necessary to rectify the cubes into that form.


Collapsing Cubes
================
Currently, the ``specutils`` library does not support loading three dimensional
spectral cubes and manipulating their non-spectral axes, so parsing a spectral
cube to a :class:`~specutils.Spectrum1D` object requires some programmatic
manipulation using the ``spectral-cube`` library. Below is an example of
loading a spectral cube object, collapsing the cube, and creating a
:class:`~specutils.Spectrum1D` object.

Loading a cube
==============

We'll use the example cube data found in the ``spectral-cube`` package, and
load the data from the repository directly into a new ``SpectralCube``
object::

    >>> from astropy.utils.data import download_file
    >>> from spectral_cube import SpectralCube
    >>> sc = SpectralCube.read(download_file("https://github.com/radio-astro-tools/spectral-cube/raw/master/spectral_cube/tests/data/example_cube.fits"), format='fits')
    >>> sc
    SpectralCube with shape=(7, 4, 3) and unit=Jy / beam:
    n_x:      3  type_x: RA---ARC  unit_x: deg    range:    52.231466 deg:   52.231544 deg
    n_y:      4  type_y: DEC--ARC  unit_y: deg    range:    31.243639 deg:   31.243739 deg
    n_s:      7  type_s: VRAD      unit_s: m / s  range:    14322.821 m / s:   14944.909 m / s

We see that the example cube is small, with a shape of ``(7, 3, 4)``; two
spatial axes ``RA`` and ``DEC``, and one spectral axis that's in velocity
units.

Extracting the spectral data
============================

Now that we have a data cube, our mission is to collapse this cube along the
two spatial axes to get the representative spectral data. We can then parse the
result into a :class:`~specutils.Spectrum1D` object. ``SpectralCube`` can take
and apply arbitrary functions over axes of the data, in this case we can apply
a simple median function over the spatial axes::

    >>> import numpy as np
    >>> sc_spec = sc.apply_numpy_function(np.mean, axis=(1,2), projection=True)
    >>> sc_spec.shape
    (7,)

The ``projection`` keyword returns a lower dimesional object that retains the WCS
information that we'll need to use when constructing the `Spectrum1D` object.

Creating a spectrum object
==========================

With our data parsed and extracted, let's go ahead and create our
:class:`~specutils.Spectrum1D` object::

    >>> from specutils import Spectrum1D
    >>> spec = Spectrum1D(flux=sc_spec.data[:] * sc_spec.unit, wcs=sc_spec.wcs)

Let's plot and see what this resulting spectrum looks like::

    >>> from matplotlib import pyplot as plt
    >>> from astropy.visualization import quantity_support
    >>> quantity_support()  # for getting units on the axes below  # doctest: +IGNORE_OUTPUT

    >>> f, ax = plt.subplots()  # doctest: +IGNORE_OUTPUT
    >>> ax.step(spec.spectral_axis, spec.flux) # doctest: +IGNORE_OUTPUT +REMOTE_DATA
