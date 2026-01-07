=============
WCS Utilities
=============

The `~specutils.utils.wcs_utils` module has functions for converting spectral values
between air and vacuum, as well as a function for calculating the refractive index
of air, which is used in the air to vacuum conversions. There are multiple methods
available for the refractive index calculation, most of which agree with each other
to a high level of precision. However, the Griesen (2006) equation seems to have a small
offset from the others of about 0.07 Angstrom (at the wavelengths shown in the plot below)
and thus in specutils 1.17.0 we have changed from using that method as the default to
using the Morton (2000) equation by default, which is consistent with the IAU standard.

The source of the discrepancy between Griesen (2006) and the other methods seems to be
that Griesen (2006) assumes an air temperature of 0C, vs an apparent assumed air temperature of 15C
for the other methods (thanks to Jon Holtzman for investigating this). Users should keep this
in mind when choosing which method to use for air to vacuum conversions.

.. image:: air_to_vac_offset.png
   :alt: Results of air to vac conversion at optical wavelengths.

The downside of all methods but Griesen (2006) is that they have mathematical singularities
in the far UV, and thus are only valid at wavelengths longer than 200 nm. The specutils
conversion functions will raise an error if these methods are used for wavelengths shorter
than this limit.

.. image:: refractive_index_singularities.png
   :alt: Demonstration of singularities in the equations for refractive index.

For additional context and discussion, see https://github.com/astropy/specutils/issues/1162.

Reference/API
-------------
.. automodapi:: specutils.utils.wcs_utils
    :no-heading:

    :skip: GWCS
    :skip: Tabular1D
    :skip: SpectralGWCS
