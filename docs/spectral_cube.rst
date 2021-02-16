###########################
Working with Spectral Cubes
###########################

Spectral cubes can be read directly with :class:`~specutils.Spectrum1D`.
A specific example of this is demonstrated here.

Note that the workflow described here is for spectral cubes that are rectified
such that one of the axes is entirely spectral and all the spaxels have the same
``spectral_axis`` values - i.e., case 2 in :ref:`specutils-representation-overview`.
For less-rectified cubes, pre-processing steps (not addressed by specutils at the
time of this writing) will be necessary to rectify the cubes into that form.


Loading a cube
==============

We'll use the example cube data ``MaNGA cube``, and load the data from the
repository directly into a new ``Spectrum1D`` object::

    >>> from astropy.utils.data import download_file
    >>> from spectral_cube import SpectralCube
    >>> from specutils.spectra import Spectrum1D
    >>> filename = "file:/Users/busko/Downloads/mangacube_cont_sub.fits"
    >>> # filename = "https://stsci.box.com/s/rhhzj5jw3pddvs6u455jhbppn0b1bw5g"
    >>> file = download_file(filename)
    >>> sc = Spectrum1D.read(file, format='MaNGA cube')
    >>> sc.shape
    (74, 74, 4563)
    >>> sc
    <Spectrum1D(flux=<Quantity [[[-1.8, -1.8, -1.8, ..., -1.8, -1.8, -1.8],
            [-1.8, -1.8, -1.8, ..., -1.8, -1.8, -1.8],
            [-1.8, -1.8, -1.8, ..., -1.8, -1.8, -1.8],
            ...,
            [-1.8, -1.8, -1.8, ..., -1.8, -1.8, -1.8],
            [-1.8, -1.8, -1.8, ..., -1.8, -1.8, -1.8],
            [-1.8, -1.8, -1.8, ..., -1.8, -1.8, -1.8]],

           [[-1.8, -1.8, -1.8, ...
           ...,
            [-1.8, -1.8, -1.8, ..., -1.8, -1.8, -1.8],
            [-1.8, -1.8, -1.8, ..., -1.8, -1.8, -1.8],
            [-1.8, -1.8, -1.8, ..., -1.8, -1.8, -1.8]]] 1e-17 erg / (Angstrom cm2 s spaxel)>, spectral_axis=<SpectralAxis [ 3621.59598486,  3622.42998417,  3623.26417553, ..., 10349.03843826,
   10351.42166679, 10353.80544415] Angstrom>, uncertainty=InverseVariance([[[0., 0., 0., ..., 0., 0., 0.],
                  [0., 0., 0., ..., 0., 0., 0.],
                  [0., 0., 0., ..., 0., 0., 0.],
                  ...,
                  [0., 0., 0., ..., 0., 0., 0.],
                  [0., 0., 0., ..., 0., 0., 0.],
                  [0., 0., 0., ..., 0., 0., 0.]]]))>


We see that the example has shape of ``(74, 74, 4563)``; two spatial axes in
a 74 X 74 configuration, plus one spectral axis that's in Angstrom units.
