import copy

import numpy as np
from astropy import units as u
from astropy.modeling.models import Shift
from astropy.modeling.tabular import Tabular1D
from gwcs import WCS as GWCS
from gwcs import coordinate_frames as cf


class SpectralGWCS(GWCS):
    """
    This is a placeholder lookup-table GWCS created when a :class:`~specutils.Spectrum1D` is
    instantiated with a ``spectral_axis`` and no WCS.
    """
    def __init__(self, *args, **kwargs):
        self.original_unit = kwargs.pop("original_unit", "")
        super().__init__(*args, **kwargs)

    def copy(self):
        """
        Return a shallow copy of the object.

        Convenience method so user doesn't have to import the
        :mod:`copy` stdlib module.

        .. warning::
            Use `deepcopy` instead of `copy` unless you know why you need a
            shallow copy.
        """
        return copy.copy(self)

    def deepcopy(self):
        """
        Return a deep copy of the object.

        Convenience method so user doesn't have to import the
        :mod:`copy` stdlib module.
        """
        return copy.deepcopy(self)

    def pixel_to_world(self, *args, **kwargs):
        if self.original_unit == '':
            return u.Quantity(super().pixel_to_world_values(*args, **kwargs))
        return super().pixel_to_world(*args, **kwargs).to(
            self.original_unit, equivalencies=u.spectral())


def refraction_index(wavelength, method='Morton2000', co2=None):
    """
    Calculates the index of refraction of dry air at standard temperature
    and pressure, at different wavelengths, using different methods.

    Parameters
    ----------
    wavelength : `Quantity` object (number or sequence)
        Vacuum wavelengths with an astropy.unit.

    method : str, optional
        Method used to convert wavelengths. Options are:

        'Morton2000' (default) - from Morton (2000, ApJS, 130, 403), eqn 8. Used by VALD,
            the Vienna Atomic Line Database. Very similar to Edlen (1966).

        'Griesen2006' - from Greisen et al. (2006, A&A 446, 747),
            eqn. 65, standard used by International Union of Geodesy and Geophysics

        'Edlen1953' - from Edlen (1953, J. Opt. Soc. Am, 43, 339). Standard
            adopted by IAU (resolution No. C15, Commission 44, XXI GA, 1991),
            which refers to Oosterhoff (1957) that uses Edlen (1953). Also used
            by Morton (1991, ApJS, 77, 119), which is frequently cited as IAU source.

        'Edlen1966' - from Edlen (1966, Metrologia 2, 71), rederived constants
            from optical and near UV data.

        'PeckReeder1972' - from Peck & Reeder (1972, J. Opt. Soc. 62), derived
            from additional infrared measurements (up to 1700 nm).

        'Ciddor1996' - from Ciddor (1996, Appl. Opt. 35, 1566). Based on
            Peck & Reeder (1972), but updated to account for the changes in
            the international temperature scale and adjust the results for
            CO2 concentration. Arguably most accurate conversion available.

        Note that all options except for 'Griesen2006' have singularities in the far
        UV. 'Griesen2006' gives values that are slightly inconsistent with the
        other methods (~0.07 Angstrom difference at visible wavelengths), but it is
        the best option in the FUV due to the mathematical singularities in the others.
        See https://specutils.readthedocs.io/en/latest/wcs_utils.html for more detail, or
        https://github.com/astropy/specutils/issues/1162 for additional context.

    co2 : number, optional
        CO2 concentration in ppm. Only used for method='Ciddor1996'. If not
        given, a default concentration of 450 ppm is used.

    Returns
    -------
    refr : number or sequence
        Index of refraction at each given air wavelength.
    """
    VALID_METHODS = ['Griesen2006', 'Edlen1953', 'Edlen1966', 'Morton2000',
                     'PeckReeder1972', 'Ciddor1996']
    assert isinstance(method, str), 'method must be a string'
    if method != 'Griesen2006' and wavelength.min() < 200 * u.nm:
        raise ValueError("The chosen method is invalid for wavelengths below 250 nm."
                         " 'Griesen2006' is the only option for this wavelength range -"
                         " see the specutils.utils.wcs_utils.refraction_index docstring"
                         " for more detail.")
    method = method.lower()
    sigma2 = (1 / wavelength.to(u.um).value)**2
    if method == 'griesen2006':
        refr = 1e-6 * (287.6155 + 1.62887 * sigma2 + 0.01360 * sigma2**2)
    elif method == 'edlen1953':
        refr = 6.4328e-5 + 2.94981e-2 / (146 - sigma2) + 2.5540e-4 / (41 - sigma2)
    elif method == 'edlen1966':
        refr = 8.34213e-5 + 2.406030e-2 / (130 - sigma2) + 1.5997e-4 / (38.9 - sigma2)
    elif method == 'morton2000':
        refr = 8.34254e-5 + 2.406147e-2 / (130 - sigma2) + 1.5998e-4 / (38.9 - sigma2)
    elif method == 'peckreeder1972':
        refr = 5.791817e-2 / (238.0185 - sigma2) + 1.67909e-3 / (57.362 - sigma2)
    elif method == 'ciddor1996':
        refr = 5.792105e-2 / (238.0185 - sigma2) + 1.67917e-3 / (57.362 - sigma2)
        if co2:
            refr *= 1 + 0.534e-6 * (co2 - 450)
    else:
        raise ValueError("Method must be one of " + ", ".join(VALID_METHODS))
    return refr + 1


def vac_to_air(wavelength, method='Morton2000', co2=None):
    """
    Converts vacuum to air wavelengths using different methods.

    Parameters
    ----------
    wavelength : `Quantity` object (number or sequence)
        Vacuum wavelengths with an astropy.unit.
    method : str, optional
        One of the methods in refraction_index(), default is 'Morton2000'.
    co2 : number, optional
        Atmospheric CO2 concentration in ppm. Only used for method='Ciddor1996'.
        If not given, a default concentration of 450 ppm is used.

    Returns
    -------
    air_wavelength : `Quantity` object (number or sequence)
        Air wavelengths with the same unit as wavelength.
    """
    refr = refraction_index(wavelength, method=method, co2=co2)
    return wavelength / refr


def air_to_vac(wavelength, scheme='inversion', method='Morton2000', co2=None,
               precision=1e-12, maxiter=30):
    """
    Converts air to vacuum wavelengths using different methods.

    Parameters
    ----------
    wavelength : `Quantity` object (number or sequence)
        Air wavelengths with an astropy.unit.

    scheme : str, optional
        How to convert from vacuum to air wavelengths. Options are:

            * 'inversion' (default) - result is simply the inversion (1 / n) of the
              refraction index of air. Griesen et al. (2006) report that the error
              in naively inverting is less than 10^-9.

            * 'Piskunov' - uses an analytical solution derived by Nikolai Piskunov
              and used by the Vienna Atomic Line Database (VALD).

            * 'iteration' - uses an iterative scheme to invert the index of refraction.

    method : str, optional
        Only used if scheme is 'inversion' or 'iteration'. One of the methods
        in `~specutils.utils.wcs_utils.refraction_index`, default is 'Morton2000'

    co2 : number, optional
        Atmospheric CO2 concentration in ppm. Only used if scheme='inversion' and
        method='Ciddor1996'. If not given, a default concentration of 450 ppm is used.

    precision : float
        Maximum fractional value in refraction conversion beyond at which iteration will
        be stopped. Only used if scheme='iteration'.

    maxiter : integer
        Maximum number of iterations to run. Only used if scheme='iteration'.

    Returns
    -------
    vac_wavelength : `Quantity` object (number or sequence)
        Vacuum wavelengths with the same unit as wavelength.
    """
    VALID_SCHEMES = ['inversion', 'iteration', 'piskunov']
    assert isinstance(scheme, str), 'scheme must be a string'
    scheme = scheme.lower()
    if scheme == 'inversion':
        refr = refraction_index(wavelength, method=method, co2=co2)
    elif scheme == 'piskunov':
        wlum = wavelength.to(u.angstrom).value
        sigma2 = (1e4 / wlum)**2
        refr = (8.336624212083e-5 + 2.408926869968e-2 / (130.1065924522 - sigma2) +
                1.599740894897e-4 / (38.92568793293 - sigma2)) + 1
    elif scheme == 'iteration':
        # Refraction index is a function of vacuum wavelengths.
        # Iterate to get index of refraction that gives air wavelength that
        # is consistent with the reverse transformation.
        counter = 0
        result = wavelength.copy()
        refr = refraction_index(wavelength, method=method, co2=co2)
        while True:
            counter += 1
            diff = wavelength * refr - result
            if abs(diff.max().value) < precision:
                break
            if counter > maxiter:
                raise RuntimeError("Reached maximum number of iterations "
                                   "without reaching desired precision level.")
            result += diff
            refr = refraction_index(result, method=method, co2=co2)
    else:
        raise ValueError("Method must be one of " + ", ".join(VALID_SCHEMES))
    return wavelength * refr


def air_to_vac_deriv(wavelength, method='Griesen2006'):
    """
    Calculates the derivative d(wave_vacuum) / d(wave_air) using different
    methods.

    Parameters
    ----------
    wavelength : `Quantity` object (number or sequence)
        Air wavelengths with an astropy.unit.

    method : str, optional
        Method used to convert wavelength derivative. Options are:
        'Griesen2006' (default) - from Greisen et al. (2006, A&A 446, 747), eqn. 66.

    Returns
    -------
    wave_deriv : `Quantity` object (number or sequence)
        Derivative d(wave_vacuum) / d(wave_air).
    """
    assert method.lower() == 'griesen2006', "Only supported method is 'Griesen2006'"
    wlum = wavelength.to(u.um).value
    return (1 + 1e-6 * (287.6155 - 1.62887 / wlum**2 - 0.04080 / wlum**4))


def gwcs_from_array(array):
    """
    Create a new WCS from provided tabular data. This defaults to being
    a GWCS object.
    """
    orig_array = u.Quantity(array)

    coord_frame = cf.CoordinateFrame(naxes=1,
                                     axes_type=('SPECTRAL',),
                                     axes_order=(0,))
    spec_frame = cf.SpectralFrame(unit=array.unit, axes_order=(0,))

    # In order for the world_to_pixel transformation to automatically convert
    # input units, the equivalencies in the look up table have to be extended
    # with spectral unit information.
    SpectralTabular1D = type("SpectralTabular1D", (Tabular1D,),
                             {'input_units_equivalencies': {'x0': u.spectral()}})

    forward_transform = SpectralTabular1D(np.arange(len(array)),
                                          lookup_table=array)
    # If our spectral axis is in descending order, we have to flip the lookup
    # table to be ascending in order for world_to_pixel to work.
    if len(array) == 0 or array[-1] > array[0]:
        forward_transform.inverse = SpectralTabular1D(
            array, lookup_table=np.arange(len(array)))
    else:
        forward_transform.inverse = SpectralTabular1D(
                array[::-1], lookup_table=np.arange(len(array))[::-1])

    tabular_gwcs = SpectralGWCS(original_unit = orig_array.unit,
                                forward_transform=forward_transform,
                                input_frame=coord_frame,
                                output_frame=spec_frame)

    # Store the intended unit from the origin input array
    #     tabular_gwcs._input_unit = orig_array.unit

    return tabular_gwcs


def gwcs_slice(self, item):
    """
    This is a bit of a hack in order to fix the slicing of the WCS
    in the spectral dispersion direction.  The NDData slices properly
    but the spectral dispersion result was not.

    There is code slightly downstream that sets the *number* of entries
    in the dispersion axis, this is just needed to shift to the correct
    starting element.

    When WCS gets the ability to do slicing then we might be able to
    remove this code.
    """
    # Create shift of x-axis
    if isinstance(item, int):
        shift = item
    elif isinstance(item, slice):
        shift = item.start
    else:
        raise TypeError('Unknown index type {}, must be int or slice.'.format(item))

    # Create copy as we need to modify this and return it.
    new_wcs = copy.deepcopy(self)

    if shift == 0:
        return new_wcs

    shifter = Shift(shift)

    # Get the current forward transform
    forward = new_wcs.forward_transform

    # Set the new transform
    new_wcs.set_transform(new_wcs.input_frame,
                          new_wcs.output_frame,
                          shifter | forward)

    return new_wcs
