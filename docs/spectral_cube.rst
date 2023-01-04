###########################
Working with Spectral Cubes
###########################

Spectral cubes can be read directly with :class:`~specutils.Spectrum1D`.
A specific example of this is demonstrated here. In addition to the functions
demonstrated below, as a subclass of `NDCube <https://github.com/sunpy/ndcube>`_,
:class:`~specutils.Spectrum1D` also inherits useful methods for e.g. cropping
based on combinations of world and spectral coordinates. Most of the
functionality inherited from `~ndcube.NDCube` requires initializing the
``Spectrum1D`` object with a WCS describing the coordinates for all axes of
the data.

Note that the workflow described here is for spectral cubes that are rectified
such that one of the axes is entirely spectral and all the spaxels have the same
``spectral_axis`` values (i.e., case 2 in :ref:`specutils-representation-overview`).
For less-rectified cubes, pre-processing steps (not addressed by specutils at the
time of this writing) will be necessary to rectify the cubes into that form.
Note, though, that such cubes can be stored in specutils data structures (cases
3 and 4 in :ref:`specutils-representation-overview`), which support *some* of
these behaviors, even though the fulkl set of tools do not yet apply.


Loading a cube
==============

We'll use a ``MaNGA cube`` for our example, and load the data from the
repository directly into a new ``Spectrum1D`` object:

.. code-block:: python

    >>> from astropy.utils.data import download_file
    >>> from specutils.spectra import Spectrum1D
    >>> filename = "https://stsci.box.com/shared/static/28a88k1qfipo4yxc4p4d40v4axtlal8y.fits"
    >>> file = download_file(filename, cache=True) # doctest: +REMOTE_DATA
    >>> sc = Spectrum1D.read(file, format='MaNGA cube') # doctest: +REMOTE_DATA


The cube has  74x74 spaxels with 4563 spectral axis points in each one:

.. code-block:: python

    >>> sc.shape # doctest: +REMOTE_DATA
    (74, 74, 4563)


Print the contents of 3 spectral axis points in a 3x3 spaxel array:

.. code-block:: python

    >>> sc[30:33,30:33,2000:2003] # doctest: +REMOTE_DATA
    <Spectrum1D(flux=<Quantity [[[0.48920232, 0.4987253 , 0.5098349 ],
                [0.493365  , 0.4964812 , 0.5223962 ],
                [0.49446177, 0.4909543 , 0.5304416 ]],
    <BLANKLINE>
               [[0.53114057, 0.53538376, 0.5467784 ],
                [0.53761804, 0.533159  , 0.554437  ],
                [0.5470889 , 0.54905874, 0.57109433]],
    <BLANKLINE>
               [[0.5599331 , 0.554316  , 0.5618426 ],
                [0.5763055 , 0.5668046 , 0.5774939 ],
                [0.59571505, 0.60118765, 0.59942234]]] 1e-17 erg / (Angstrom cm2 s spaxel)>, spectral_axis=<SpectralAxis
       (observer to target:
          radial_velocity=0.0 km / s
          redshift=0.0)
      [5.73984286e-07, 5.74116466e-07, 5.74248676e-07] m>, uncertainty=InverseVariance([[[4324.235 , 4326.87  , 4268.985 ],
                      [5128.3867, 5142.5005, 4998.457 ],
                      [4529.9097, 4545.8345, 4255.305 ]],
    <BLANKLINE>
                     [[4786.163 , 4811.216 , 4735.3135],
                      [4992.71  , 5082.1294, 4927.881 ],
                      [4992.9683, 5046.971 , 4798.005 ]],
    <BLANKLINE>
                     [[4831.2236, 4887.096 , 4806.84  ],
                      [3895.8677, 4027.9104, 3896.0195],
                      [4521.258 , 4630.997 , 4503.0396]]]))>


Spectral slab extraction
========================

The `~specutils.manipulation.spectral_slab` function can be used to extract
spectral regions from the cube.

.. code-block:: python

    >>> import astropy.units as u
    >>> from specutils.manipulation import spectral_slab
    >>> ss = spectral_slab(sc, 5000.*u.AA, 5003.*u.AA) # doctest: +REMOTE_DATA
    >>> ss.shape  # doctest: +REMOTE_DATA
    (74, 74, 3)
    >>> ss[30:33,30:33,::] # doctest: +REMOTE_DATA
    <Spectrum1D(flux=<Quantity [[[0.6103081 , 0.95697385, 1.0791174 ],
                [0.5663384 , 0.8872061 , 1.0814004 ],
                [0.520966  , 0.7819859 , 1.024845  ]],
    <BLANKLINE>
               [[0.64514536, 0.96376216, 1.083235  ],
                [0.6112465 , 0.89025146, 1.058679  ],
                [0.56316894, 0.77895504, 0.99165994]],
    <BLANKLINE>
               [[0.65954393, 0.9084677 , 0.9965009 ],
                [0.6255246 , 0.84401435, 0.9930112 ],
                [0.59066033, 0.762025  , 0.9361185 ]]] 1e-17 erg / (Angstrom cm2 s spaxel)>, spectral_axis=<SpectralAxis
       (observer to target:
          radial_velocity=0.0 km / s
          redshift=0.0)
      [5.00034537e-07, 5.00149688e-07, 5.00264865e-07] m>, uncertainty=InverseVariance([[[3449.242 , 2389.292 , 2225.105 ],
                      [4098.7485, 2965.88  , 2632.497 ],
                      [3589.92  , 2902.7622, 2292.3823]],
    <BLANKLINE>
                     [[3563.3342, 2586.58  , 2416.039 ],
                      [4090.8855, 3179.1702, 2851.823 ],
                      [4158.919 , 3457.0115, 2841.1965]],
    <BLANKLINE>
                     [[3684.6013, 3056.2   , 2880.6592],
                      [3221.7888, 2801.3518, 2525.541 ],
                      [3936.68  , 3461.534 , 3047.6135]]]))>


Spectral Bounding Region
========================

The `~specutils.manipulation.extract_bounding_spectral_region` function can be used to
extract the bounding region that encompases a set of disjoint `~specutils.SpectralRegion`
instances, or a composite instance of `~specutils.SpectralRegion` that contains
disjoint sub-regions.

.. code-block:: python

    >>> from specutils import SpectralRegion
    >>> from specutils.manipulation import extract_bounding_spectral_region
    >>> composite_region = SpectralRegion([(5000*u.AA, 5002*u.AA), (5006*u.AA, 5008.*u.AA)])
    >>> sub_spectrum = extract_bounding_spectral_region(sc, composite_region) # doctest: +REMOTE_DATA
    >>> sub_spectrum.spectral_axis  # doctest: +REMOTE_DATA +FLOAT_CMP
    <SpectralAxis
       (observer to target:
          radial_velocity=0.0 km / s
          redshift=0.0)
      [5.00034537e-07, 5.00149688e-07, 5.00264865e-07, 5.00380068e-07,
       5.00495298e-07, 5.00610555e-07, 5.00725838e-07] m>


Moments
=======

The `~specutils.analysis.moment` function can be used to compute moments of any order
along one of the cube's axes. By default, ``axis=-1``, which computes moments
along the spectral axis (remember that the spectral axis is always last in a
:class:`~specutils.Spectrum1D`).

.. code-block:: python

    >>> from specutils.analysis import moment
    >>> m = moment(sc, order=1) # doctest: +REMOTE_DATA
    >>> m.shape  # doctest: +REMOTE_DATA
    (74, 74)
    >>> m[30:33,30:33]  # doctest: +REMOTE_DATA +FLOAT_CMP
    <Quantity [[6.45261317e-07, 6.46265069e-07, 6.48128166e-07],
               [6.46467930e-07, 6.47941283e-07, 6.51460998e-07],
               [6.48672775e-07, 6.52631872e-07, 6.56733087e-07]] m>

Use Case
========

Example of computing moment maps for specific wavelength ranges in a
cube, using `~specutils.manipulation.spectral_slab` and
`~specutils.analysis.moment`.

.. plot::
    :include-source:
    :align: center
    :context: close-figs

    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.utils.data import download_file
    from specutils import Spectrum1D, SpectralRegion
    from specutils.analysis import moment
    from specutils.manipulation import spectral_slab

    filename = "https://stsci.box.com/shared/static/28a88k1qfipo4yxc4p4d40v4axtlal8y.fits"
    fn = download_file(filename, cache=True)
    spec1d = Spectrum1D.read(fn)

    # Extract H-alpha sub-cube for moment maps using spectral_slab
    subspec = spectral_slab(spec1d, 6745.*u.AA, 6765*u.AA)
    ha_wave = subspec.spectral_axis

    # Extract wider sub-cube covering H-alpha and [N II] using spectral_slab
    subspec_wide = spectral_slab(spec1d, 6705.*u.AA, 6805*u.AA)
    ha_wave_wide= subspec_wide.spectral_axis

    # Convert flux density to microJy and correct negative flux offset for
    # this particular dataset
    ha_flux = (np.sum(subspec.flux.value, axis=(1,2)) + 0.0093) * 1.0E-6*u.Jy
    ha_flux_wide = (np.sum(subspec_wide.flux.value, axis=(1,2)) + 0.0093) * 1.0E-6*u.Jy

    # Compute moment maps for H-alpha line
    moment0_halpha = moment(subspec, order=0)
    moment1_halpha = moment(subspec, order=1)

    # Convert moment1 from AA to velocity
    # H-alpha is redshifted to 6755 AA for this galaxy
    print(moment1_halpha[40,40])
    vel_map = 3.0E5 * (moment1_halpha.value - 6.755E-7) / 6.755E-7

    # Plot results in 3 panels (subspec_wide,  H-alpha line flux, H-alpha velocity map)
    f,(ax1,ax2,ax3) = plt.subplots(1, 3, figsize=(15, 5))
    ax1.plot(ha_wave_wide, (ha_flux_wide)*1000.)
    ax1.set_xlabel('Angstrom', fontsize=14)
    ax1.set_ylabel('uJy', fontsize=14)
    ax1.tick_params(axis="both", which='major', labelsize=14, length=8, width=2, direction='in', top=True, right=True)
    ax2.imshow(moment0_halpha.value, origin='lower')
    ax2.set_title('moment = 0')
    ax2.set_xlabel('x pixels', fontsize=14)
    ax3.imshow(vel_map, vmin=-100., vmax=100., cmap='rainbow', origin='lower')
    ax3.set_title('moment = 1')
    ax3.set_xlabel('x pixels', fontsize=14)
