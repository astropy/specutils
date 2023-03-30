1.10.0 (unreleased)
-------------------

New Features
^^^^^^^^^^^^

Bug Fixes
^^^^^^^^^

- Fix bug in ``fit_lines`` which gave unexpected outputs from the ``get_fit_info``
  and ``ignore_units`` keyword arguments [#1030]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- JWST s2d and s3d readers now requires the optional dependency, ``stdatamodels``,
  which user has to install separately. [#1038]

- ASDF tag for Spectrum1D is now compatible with ASDF v3.
  As a result, minversion of ``asdf`` has been bumped to 2.14.
  Redundant ASDF schema for ``SpectralCoord`` is removed.
  It also now supports ``mask`` serialization. [#1042, #1053]

- JWST X1D reader will no longer raise a ``UnitWarning`` for surface brightness
  error. [#1050]

1.9.1 (2022-11-22)
------------------

Bug Fixes
^^^^^^^^^

- Add and subtract operations on Spectrum1D now allow for other operand's class
  to handle the arithmetic if that class has special handling. [#988]

1.9.0 (2022-10-18)
------------------

Bug Fixes
^^^^^^^^^

- Fix bug in fitting with weights if weights argument is set to 'unc'. [#979]
- Fix bug in JWST reader which caused multi-extension files to load only the
  primary HDU [#982]
- Implemented conversion to expected uncertainty type in a few functions that
  were still just assuming the uncertainty was the correct type. [#984]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Bumped astropy minimum version to 5.1. [#984]

1.8.1 (2022-09-09)
------------------

Bug Fixes
^^^^^^^^^

- Arithmetic with constants and Spectrum1D now works in either order. [#964]

- Fixed uncertainty propagation in FluxConservingResampler. [#976]

1.8.0 (2022-08-22)
------------------

New Features
^^^^^^^^^^^^

- Implemented uncertainty propagation for analysis functions. [#938, #939, #961, #968]

- Model fitting with ``fit_lines`` now returns uncertainties from the underlying scipy
  fitter by default. [#962]

Bug Fixes
^^^^^^^^^

- Fixed a bug with moment map orders greater than 1 not being able to handle
  cubes with non-square spatial dimensions. [#970]

- Added a workaround for reading JWST IFUs with incorrect GWCS. [#973]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- The Spectrum1D redshift and radial_velocity attribute setters were deprecated
  in favor of the more explicit set_redshift_to, shift_spectrum_to, and
  set_radial_velocity_to methods. [#946, #943]

- ``estimate_line_parameters`` now calculates estimates based on the selected
  region, rather than the entire spectrum. [#962]

1.7.0 (2022-02-21)
------------------

New Features
^^^^^^^^^^^^

Bug Fixes
^^^^^^^^^

- Fixed ``spectral_slab`` crashing when ``spectral_axis`` has unit of pixels and
  the bounds are also defined in the unit of pixels. [#926]

- Fixed resulting ``spectral_axis`` containing NaN when a cube is passed into
  ``Spectrum1D`` without WCS nor spectral axis and the spatial-spatial dimension
  is smaller than spectral dimension. [#926]

- Fixed WCS not accurately reflecting the updated spectral axis after slicing a 
  ``Spectrum1D``. [#918]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Logger usage is removed. Warnings now issued using Python ``warnings`` module.
  This enables more granular warning control for downstream packages. [#922]

1.6.0 (2022-01-27)
------------------

New Features
^^^^^^^^^^^^

- Add collapse methods to Spectrum1D. [#904, #906]

- SpectralRegion and Spectrum1D now allow descending (in wavelength space) as
  well as ascending spectral values. [#911]

1.5.0
-----

New Features
^^^^^^^^^^^^

- Convolution-based smoothing will now apply a 1D kernel to multi-dimensional fluxes
  by convolving along the spectral axis only, rather than raising an error. [#885]

- ``template_comparison`` now handles ``astropy.nddata.Variance`` and 
  ``astropy.nddata.InverseVariance`` uncertainties instead of assuming 
  the uncertainty is standard deviation. [#899]

Bug Fixes
^^^^^^^^^

- Speed up JWST s3d loader and reduce memory usage. [#874]

- ``SpectralRegion`` can now handle pixels. [#886]

- Fix bug where ``template_comparison`` would return the wrong chi2 value. [#872]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``fit_lines`` now makes use of unit support in ``astropy.modeling``. [#891]

- ``Spectrum1D.with_spectral_units`` now attempts to fall back on the ``spectral_axis``
  units if units could not be retrieved from the WCS. [#892]

- ``ndcube`` package pin updated to released version (2.0). [#897]

- Minor changes for astropy 5.0 compatibility. [#895]

1.4.1
-----

Bug Fixes
^^^^^^^^^

- Fix JWST s3d loader. [#866]

1.4
---

New Features
^^^^^^^^^^^^

- Allow overriding existing astropy registry elements. [#861]

- ``Spectrum1D`` will now swap the spectral axis with the last axis on initialization 
  if it can be identified from the WCS and is not last, rather than erroring. [#654, #822]

Bug Fixes
^^^^^^^^^

- Change loader priorities so survey loaders always override generic ones. [#860]

- Handle "FLUX_ERROR" header keyword in addition to "ERROR" in JWST loader. [#856]


Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``Spectrum1D`` now subclasses ``NDCube`` instead of ``NDDataRef``. [#754, #822, #850]

1.3.1
-----

New Features
^^^^^^^^^^^^

- Add ``SpectrumList`` loader for set of JWST _x1d files. [#838]

Bug Fixes
^^^^^^^^^

- Handle new ``astropy.units.PhysicalType`` class added in astropy 4.3. [#833]
- Handle case of WCS with None values in ``world_axis_physical_types`` when
  initializing Spectrum1D. [#839]
- Fix bug in apStar loader. [#839]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Improve continuum flux calculation in ``equivalent_width``. [#843]

1.3
---

New Features
^^^^^^^^^^^^

- Added ability to slice ``Spectrum1D`` with spectral axis values. [#790]

- Added ability to replace a section of a spectrum with a spline or model fit. [#782]

Bug Fixes
^^^^^^^^^

- Fix infinite recursion when unpickling a ``QuantityModel``. [#823]

- Changed positional to keyword arguments in ``fit_continuum``. [#806]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Fix inaccuracy about custom loading in docs. [#819]

- Use non-root logger to prevent duplicate messages. [#810]

- Removed unused astropy config code. [#805]

1.2
---

New Features
^^^^^^^^^^^^

- Add support for reading IRAF MULTISPEC format with non-linear 2D WCS into
  ``SpectrumCollection`` to default_loaders. [#708]

- ``SpectralRegion`` objects can now be created from the ``QTable``
  object returned from the line finding rountines. [#759]

- Include new 6dFGS loaders. [#734]

- Include new OzDES loaders. [#764]

- Include new GAMA survey loaders. [#765]

- Include new GALAH loaders. [#766]

- Include new WiggleZ loaders. [#767]

- Include new 2dF/AAOmega loaders. [#768]

- Add loader to handle IRAF MULTISPEC non-linear 2D WCS. [#708]

- Add ability to extract minimum bounding regions of ``SpectralRegion`` objects. [#755]

- Implement new moment analysis function for specutils objects. [#758]

- Add new spectral slab extraction functionality. [#753]

- Include new loaders for AAT and other Australian surveys. [#719]

- Improve docstrings and intialization of ``SpectralRegion`` objects. [#770]


Bug Fixes
^^^^^^^^^

- Fix ``extract_region`` behavior and slicing for ``Spectrum1D`` objects
  that have multi-dimensional flux arrays. Extracting a region that extends
  beyond the limits of the data no longer drops the last data point in the
  returned spectrum. [#724]

- Fixes to the jwst loaders. [#759]

- Fix handling of ``SpectralCollection`` objects in moment calculations. [#781]

- Fix issue with non-loadable x1d files. [#775]

- Fix WCS handling in SDSS loaders. [#738]

- Fix the property setters for radial velocity and redshift. [#722]

- Fix line test errors and include python 3.9 in tests. [#751]

- Fix smoothing functionality dropping spectrum meta information. [#732]

- Fix region extraction for ``Spectrum1D`` objects with multi-dimensional fluxes. [#724]

Documentation
^^^^^^^^^^^^^

- Update SDSS spectrum documentation examples. [#778]

- Include new documentation on working with ``SpectralCube`` objects. [#726, #784]

- Add documentation on spectral cube related functionality. [#783]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Improved error messages when creating ``SpectralRegion`` objects. [#759]

- Update documentation favicons and ensure color consistency. [#780]

- Remove fallback ``SpectralCoord`` code and rely on upstream. [#786]

- Move remaining loaders to use utility functions for parsing files. [#718]

- Remove unnecessary data reshaping in tabular fits writer. [#730]

- Remove astropy helpers and CI helpers dependencies. [#562]

1.1
---

New Features
^^^^^^^^^^^^

- Added writer to ``wcs1d-fits`` and support for multi-D flux arrays with
  1D WCS (identical ``spectral_axis`` scale). [#632]

- Implement ``SpectralCoord`` for ``SpectrumCollection`` objects. [#619]

- Default loaders work with fits file-like objects. [#637]

- Implement bin edge support on ``SpectralCoord`` objects using
  ``SpectralAxis`` subclass. [#645]

- Implement new 6dFGS loader. [#608]

- Implement uncertainty handling for ``line_flux``. [#669]

- Implement new 2SLAQ-LRG loader. [#633]

- Implement new 2dFGRS loader. [#695]

- Default loaders now include WCS 1D (with multi-dimensional flux handling) writer. [#632]

- Allow continuum fitting over multiple windows. [#698]

- Have NaN-masked inputs automatically update the ``mask`` appropriately. [#699]

Bug Fixes
^^^^^^^^^

- Fixed ``tabular-fits`` handling of 1D+2D spectra without WCS;
  identification and parsing of metadata and units for ``apogee``
  and ``muscles`` improved; enabled loading from file-like objects. [#573]

- Fix ASDF handling of ``SpectralCoord``. [#642]

- Preserve flux unit in ``resample1d`` for older versions of numpy. [#649]

- Fix setting the doppler values on ``SpectralCoord`` instances. [#657]

- Properly handle malformed distances in ``SkyCoord`` instances. [#663]

- Restrict spectral equivalencies to contexts where it is required. [#573]

- Fix ``from_center`` descending spectral axis handling. [#656]

- Fix factor of two error in ``from_center`` method of ``SpectralRegion`` object. [#710]

- Fix handling of multi-dimensional mask slicing. [#704]

- Fix identifier for JWST 1D loader. [#715]

Documentation
^^^^^^^^^^^^^

- Display supported loaders in specutils documentation. [#675]

- Clarify inter-relation of specutils objects in relevant docstrings. [#654]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Remove pytest runtime dependency. [#603]

- Change implementation of ``.quantity`` to ``.view`` in ``SpectralCoord``. [#614]

- Ensure underlying references point to ``SpectralCoord`` object. [#640]

- Deprecate ``spectral_axis_unit`` property. [#618]

- Backport ``SpectralCoord`` from astropy core for versions <4.1. [#674]

- Improve SDSS loaders and improve handling of extensions. [#667]

- Remove spectral cube testing utilities. [#683]

- Change local specutils directory creation behavior. [#691]

- Ensure existing manipulation and analysis functions use ``mask`` attribute. [#670]

- Improve mask handling in analysis functions. [#701]

1.0
---

New Features
^^^^^^^^^^^^

- Implement ``SpectralCoord`` object. [#524]

- Implement cross-correlation for finding redshift/radial velocity. [#544]

- Improve FITS file identification in default_loaders. [#545]

- Support ``len()`` for ``SpectrumCollection`` objects. [#575]

- Improved 1D JWST loader and allow parsing into an ``SpectrumCollection`` object. [#579]

- Implemented 2D and 3D data loaders for JWST products. [#595]

- Include documentation on how to use dust_extinction in specutils. [#594]

- Include example of spectrum shifting in docs. [#600]

- Add new default excise_regions exciser function and improve subregion handling. [#609]

- Implement use of ``SpectralCoord`` in ``Spectrum1D`` objects. [#610]

Bug Fixes
^^^^^^^^^

- Fix stacking and unit treatment in ``SpectrumCollection.from_spectra``. [#578]

- Fix spectral axis unit retrieval. [#581]

- Fix bug in subspectrum fitting. [#586]

- Fix uncertainty to weight conversion to match astropy assumptions. [#594]

- Fix warnings and log messages from ASDF serialization tests. [#597]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Remove spectral_resolution stub from Spectrum1D. [#606]


0.7
---

New Features
^^^^^^^^^^^^

- Make specutils compatible with Astropy 4.0 (breaking change). [#462]

- Remove all wcs adapter code and rely on APE14 implementation. [#462]

Bug Fixes
^^^^^^^^^

- Address ``MexicanHat1D`` name change in documentation. [#564]


0.6.1
-----

API Changes
^^^^^^^^^^^

- Resamplers now include ``extrapolation_treatment`` argument. [#537]

- Template fitting now returns an array of chi squared values for each template. [#551]

New Features
^^^^^^^^^^^^

- Masks now supported in fitting operations. [#519]

- Resamplers now support resamping beyond the edge of a spectrum using. [#537]

- New template fitting for redshift finding. [#527]

- New continuum checker to discern whether continuum is normalized or subtracted. [#538]

- Include documentation on how to achieve splicing using specutils. [#536]

- Include function to calculate masks based on S/N thresholding. [#509]

Bug Fixes
^^^^^^^^^

- Include new regions regression tests. [#345]

- Fix fitting documentation code block test. [#478]

- Fix Apogee loader to incorporate spectral axis keyword argument. [#560]

- Fix tabular fits writer and include new regression test. [#539]

- Fix dispersion attribute bug in ``Spectrum1D`` objects. [#530]

- Correctly label regression tests that require remote data. [#525]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Switch to using ``gaussian_sigma_width`` for ``Gaussian1D`` fitting estimator. [#434]

- Update documentation side bar to include page listing. [#556]

- New documentation on ``spectrum_mixin``. [#532]

- Model names are now preserved in the ``fit_lines`` operation. [#526]

- Clearer error messages for incompatible units in line fitting. [#520]

- Include travis stages in test matrix. [#515]


0.6
---

New Features
^^^^^^^^^^^^

- New redshift and radial velocity storage on `Spectrum1D` object.

- Spectral template matching including resampling.

- Error propagation in convolution smoothing.

- Sub-pixel precision for fwhm calculations.

- New spectral resampling functions.

- New IRAF data loaders.

- New FWZI calculation.

Bug Fixes
^^^^^^^^^

- Stricter intiailizer for ``Spectrum1D``.

- Correct handling of weights in line fitting.

- Array size checking in `Spectrum1D` objects.

- Fix for continuum fitting on pixel-axis dispersions.

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

0.5.3
-----

Bug Fixes
^^^^^^^^^

- Fix comparison of FITSWCS objects in arithmetic operations.

- Fix example documentation when run in python interpreter.


0.5.2 (2019-02-06)
----------------

Bug Fixes
^^^^^^^^^
- Bugfixes for astropy helpers, pep8 syntax checking, and plotting in docs [#416,#417,#419]

- All automatically generated ``SpectrumList`` loaders now have identifiers. [#440]

- ``SpectralRange.from_center`` parameters corrected after change to SpectralRange interface. [#433]

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Improve explanation on creating spectrum continua. [#420]

- Wrap IO identifier functions to ensure they always return True or False and log any errors. [#404]


0.5.1 (2018-11-29)
------------------

Bug Fixes
^^^^^^^^^

- Fixed a bug in using spectral regions that have been inverted. [#403]

- Use the pytest-remotedata plugin to control tests that require access to remote data. [#401,#408]


0.5 (2018-11-21)
----------------

This was the first release of specutils executing the
[APE14](https://github.com/astropy/astropy-APEs/blob/main/APE14.rst)
plan (i.e. the "new" specutils) and therefore intended for broad use.
