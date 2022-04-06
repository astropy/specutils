import warnings

import astropy.units as u
from astropy import constants
from astropy.utils.decorators import lazyproperty
from astropy.coordinates import SpectralCoord
import numpy as np

__all__ = ['SpectralAxis']

# We don't want to run doctests in the docstrings we inherit from Quantity
__doctest_skip__ = ['SpectralAxis.*']


class SpectralAxis(SpectralCoord):
    """
    Coordinate object representing spectral values corresponding to a specific
    spectrum. Overloads SpectralCoord with additional information (currently
    only bin edges).

    Parameters
    ----------
    bin_specification: str, optional
        Must be "edges" or "centers". Determines whether specified axis values
        are interpreted as bin edges or bin centers. Defaults to "centers".
    """

    _equivalent_unit = SpectralCoord._equivalent_unit + (u.pixel,)

    def __new__(cls, value, *args, bin_specification="centers", **kwargs):

        # Convert to bin centers if bin edges were given, since SpectralCoord
        # only accepts centers
        if bin_specification == "edges":
            bin_edges = value
            value = SpectralAxis._centers_from_edges(value)

        obj = super().__new__(cls, value, *args, **kwargs)

        if bin_specification == "edges":
            obj._bin_edges = bin_edges

        return obj

    @staticmethod
    def _edges_from_centers(centers, unit):
        """
        Calculates interior bin edges based on the average of each pair of
        centers, with the two outer edges based on extrapolated centers added
        to the beginning and end of the spectral axis.
        """
        a = np.insert(centers, 0, 2*centers[0] - centers[1])
        b = np.append(centers, 2*centers[-1] - centers[-2])
        edges = (a + b) / 2
        return edges*unit

    @staticmethod
    def _centers_from_edges(edges):
        """
        Calculates the bin centers as the average of each pair of edges
        """
        return (edges[1:] + edges[:-1]) / 2

    @lazyproperty
    def bin_edges(self):
        """
        Calculates bin edges if the spectral axis was created with centers
        specified.
        """
        if hasattr(self, '_bin_edges'):
            return self._bin_edges
        else:
            return self._edges_from_centers(self.value, self.unit)

    def with_observer_stationary_relative_to(self, frame,
                                             velocity=None,
                                             preserve_observer_frame=False):
        if self.unit is u.pixel:
            raise u.UnitsError("Cannot transform spectral coordinates in pixel units")
        super().with_observer_stationary_relative_to(frame,
                                                     velocity=velocity,
                                                     preserve_observer_frame=preserve_observer_frame)

    def with_radial_velocity_shift(self, target_shift=None, observer_shift=None):
        if self.unit is u.pixel:
            raise u.UnitsError("Cannot transform spectral coordinates in pixel units")
        return super().with_radial_velocity_shift(target_shift=target_shift,
                                                  observer_shift=observer_shift)


# ----- Copied from spectral-cube -----

DOPPLER_CONVENTIONS = {'radio': u.doppler_radio,
                       'optical': u.doppler_optical,
                       'relativistic': u.doppler_relativistic}

# These are the only linear transformations allowed
LINEAR_CTYPES = {u.doppler_optical: 'VOPT', u.doppler_radio: 'VRAD',
                 u.doppler_relativistic: 'VELO'}
LINEAR_CTYPE_CHARS = {u.doppler_optical: 'W', u.doppler_radio: 'F',
                      u.doppler_relativistic: 'V'}

ALL_CTYPES = {'speed': LINEAR_CTYPES,
              'frequency': 'FREQ',
              'length': 'WAVE'}

CTYPE_TO_PHYSICALTYPE = {'WAVE': 'length',
                         'AIR': 'air wavelength',
                         'AWAV': 'air wavelength',
                         'FREQ': 'frequency',
                         'VELO': 'speed',
                         'VRAD': 'speed',
                         'VOPT': 'speed'}

CTYPE_CHAR_TO_PHYSICALTYPE = {'W': 'length',
                              'A': 'air wavelength',
                              'F': 'frequency',
                              'V': 'speed'}
CTYPE_TO_PHYSICALTYPE.update(CTYPE_CHAR_TO_PHYSICALTYPE)

PHYSICAL_TYPE_TO_CTYPE = dict([(v,k) for k,v in
                               CTYPE_CHAR_TO_PHYSICALTYPE.items()])
PHYSICAL_TYPE_TO_CHAR = {'speed': 'V',
                         'frequency': 'F',
                         'length': 'W'}

# Used to indicate the intial / final sampling system
WCS_UNIT_DICT = {'F': u.Hz, 'W': u.m, 'V': u.m/u.s}
PHYS_UNIT_DICT = {'length': u.m, 'frequency': u.Hz, 'speed': u.m/u.s}

LINEAR_CUNIT_DICT = {'VRAD': u.Hz, 'VOPT': u.m, 'FREQ': u.Hz, 'WAVE': u.m,
                     'VELO': u.m/u.s, 'AWAV': u.m}
LINEAR_CUNIT_DICT.update(WCS_UNIT_DICT)


def _parse_velocity_convention(vc):
    if vc in (u.doppler_radio, 'radio', 'RADIO', 'VRAD', 'F', 'FREQ'):
        return u.doppler_radio
    elif vc in (u.doppler_optical, 'optical', 'OPTICAL', 'VOPT', 'W', 'WAVE'):
        return u.doppler_optical
    elif vc in (u.doppler_relativistic, 'relativistic', 'RELATIVE', 'VREL',
                'speed', 'V', 'VELO'):
        return u.doppler_relativistic


def parse_phys_type(unit):
    '''
    As of astropy 4.3, the physical type of a speed is now "speed/velocity".
    This is to parse those types and return "speed" that works with our dictionary defintions,
    and will also continue to work with previous astropy versions.
    '''
    return 'speed' if 'speed' in str(unit.physical_type) else str(unit.physical_type)


def get_rest_value_from_wcs(mywcs):
    if mywcs.wcs.restfrq:
        ref_value = mywcs.wcs.restfrq * u.Hz
        return ref_value
    elif mywcs.wcs.restwav:
        ref_value = mywcs.wcs.restwav * u.m
        return ref_value


def cdelt_derivative(crval, cdelt, intype, outtype, linear=False, rest=None):
    if intype == outtype:
        return cdelt
    elif set((outtype,intype)) == set(('length','frequency')):
        # Symmetric equations!
        return (-constants.c / crval**2 * cdelt).to(PHYS_UNIT_DICT[outtype])
    elif outtype in ('frequency','length') and 'speed' in intype:
        if linear:
            numer = cdelt * rest.to(PHYS_UNIT_DICT[outtype], u.spectral())
            denom = constants.c
        else:
            numer = cdelt * constants.c * rest.to(PHYS_UNIT_DICT[outtype], u.spectral())
            denom = (constants.c + crval)*(constants.c**2 - crval**2)**0.5
        if outtype == 'frequency':
            return (-numer/denom).to(PHYS_UNIT_DICT[outtype], u.spectral())
        else:
            return (numer/denom).to(PHYS_UNIT_DICT[outtype], u.spectral())
    elif 'speed' in outtype and intype in ('frequency','length'):

        if linear:
            numer = cdelt * constants.c
            denom = rest.to(PHYS_UNIT_DICT[intype], u.spectral())
        else:
            numer = 4 * constants.c * crval * rest.to(crval.unit, u.spectral())**2 * cdelt
            denom = (crval**2 + rest.to(crval.unit, u.spectral())**2)**2
        if intype == 'frequency':
            return (-numer/denom).to(PHYS_UNIT_DICT[outtype], u.spectral())
        else:
            return (numer/denom).to(PHYS_UNIT_DICT[outtype], u.spectral())
    elif intype == 'air wavelength':
        raise TypeError("Air wavelength should be converted to vacuum earlier.")
    elif outtype == 'air wavelength':
        raise TypeError("Conversion to air wavelength not supported.")
    else:
        raise ValueError("Invalid in/out frames")


def air_to_vac_deriv(wavelength):
    """
    Eqn 66
    """
    wlum = wavelength.to(u.um).value
    return (1 + 1e-6 * (287.6155 - 1.62887 / wlum**2 - 0.04080 / wlum**4))


def air_to_vac(wavelength):
    """
    Implements the air to vacuum wavelength conversion described in eqn 65 of
    Griesen 2006
    """
    wlum = wavelength.to(u.um).value
    return (1 + 1e-6 * (287.6155 + 1.62887 / wlum**2 + 0.01360/ wlum**4)) * wavelength


def ctype_from_vconv(ctype, unit, velocity_convention=None):
    """
    Given a CTYPE describing the current WCS and an output unit and velocity
    convention, determine the appropriate output CTYPE

    Examples
    --------
    >>> determine_ctype_from_vconv('VELO-F2V', u.Hz)
    'FREQ'
    >>> determine_ctype_from_vconv('VELO-F2V', u.m)
    'WAVE-F2W'
    >>> determine_ctype_from_vconv('FREQ', u.m/u.s)  # doctest: +SKIP
    ...
    ValueError: A velocity convention must be specified
    >>> determine_ctype_from_vconv('FREQ', u.m/u.s, velocity_convention=u.doppler_radio)
    'VRAD'
    >>> determine_ctype_from_vconv('FREQ', u.m/u.s, velocity_convention=u.doppler_optical)
    'VOPT-F2W'
    >>> determine_ctype_from_vconv('FREQ', u.m/u.s, velocity_convention=u.doppler_relativistic)
    'VELO-F2V'
    """
    unit = u.Unit(unit)

    if len(ctype) > 4:
        in_physchar = ctype[5]
    else:
        lin_cunit = LINEAR_CUNIT_DICT[ctype]

        in_physchar = PHYSICAL_TYPE_TO_CHAR[parse_phys_type(lin_cunit)]

    if parse_phys_type(unit) == 'speed':
        if velocity_convention is None and ctype[0] == 'V':
            # Special case: velocity <-> velocity doesn't care about convention
            return ctype
        elif velocity_convention is None:
            raise ValueError('A velocity convention must be specified')
        vcin = _parse_velocity_convention(ctype[:4])
        vcout = _parse_velocity_convention(velocity_convention)
        if vcin == vcout:
            return LINEAR_CTYPES[vcout]
        else:
            return "{type}-{s1}2{s2}".format(type=LINEAR_CTYPES[vcout],
                                             s1=in_physchar,
                                             s2=LINEAR_CTYPE_CHARS[vcout])

    else:
        in_phystype = CTYPE_TO_PHYSICALTYPE[in_physchar]
        if in_phystype == parse_phys_type(unit):
            # Linear case
            return ALL_CTYPES[in_phystype]
        else:
            # Nonlinear case
            out_physchar = PHYSICAL_TYPE_TO_CTYPE[parse_phys_type(unit)]
            return "{type}-{s1}2{s2}".format(type=ALL_CTYPES[parse_phys_type(unit)],
                                             s1=in_physchar,
                                             s2=out_physchar)


def convert_spectral_axis(mywcs, outunit, out_ctype, rest_value=None):
    """
    Convert a spectral axis from its unit to a specified out unit with a given output
    CTYPE.

    Only VACUUM units are supported (not air).

    Process:

    1. Convert the input unit to its equivalent linear unit.
    2. Convert the input linear unit to the output linear unit.
    3. Convert the output linear unit to the output unit.

    """
    # If the WCS includes a rest frequency/wavelength, convert it to frequency
    # or wavelength first.  This allows the possibility of changing the rest
    # frequency
    wcs_rv = get_rest_value_from_wcs(mywcs)
    inunit = u.Unit(mywcs.wcs.cunit[mywcs.wcs.spec])
    outunit = u.Unit(outunit)

    # If wcs_rv is set and speed -> speed, then we're changing the reference
    # location and we need to convert to meters or Hz first
    if ((parse_phys_type(inunit) == 'speed' and
         parse_phys_type(outunit) == 'speed' and
         wcs_rv is not None)):
        mywcs = convert_spectral_axis(mywcs, wcs_rv.unit,
                                      ALL_CTYPES[parse_phys_type(wcs_rv.unit)],
                                      rest_value=wcs_rv)
        inunit = u.Unit(mywcs.wcs.cunit[mywcs.wcs.spec])
    elif (parse_phys_type(inunit) == 'speed' and parse_phys_type(outunit) == 'speed'
          and wcs_rv is None):
        # If there is no reference change, we want an identical WCS, since
        # WCS doesn't know about units *at all*
        newwcs = mywcs.deepcopy()
        return newwcs

    in_spec_ctype = mywcs.wcs.ctype[mywcs.wcs.spec]

    # Check whether we need to convert the rest value first
    ref_value = None
    if 'speed' in parse_phys_type(outunit):
        if rest_value is None:
            rest_value = wcs_rv
            if rest_value is None:
                raise ValueError("If converting from wavelength/frequency to speed, "
                                 "a reference wavelength/frequency is required.")
        ref_value = rest_value.to(u.Hz, u.spectral())
    elif 'speed' in parse_phys_type(inunit):
        # The rest frequency and wavelength should be equivalent
        if rest_value is not None:
            ref_value = rest_value
        elif wcs_rv is not None:
            ref_value = wcs_rv
        else:
            raise ValueError("If converting from speed to wavelength/frequency, "
                             "a reference wavelength/frequency is required.")

    # If the input unit is not linearly sampled, its linear equivalent will be
    # the 8th character in the ctype, and the linearly-sampled ctype will be
    # the 6th character
    # e.g.: VOPT-F2V
    lin_ctype = (in_spec_ctype[7] if len(in_spec_ctype) > 4 else in_spec_ctype[:4])
    lin_cunit = (LINEAR_CUNIT_DICT[lin_ctype] if lin_ctype in LINEAR_CUNIT_DICT
                 else mywcs.wcs.cunit[mywcs.wcs.spec])
    in_vcequiv = _parse_velocity_convention(in_spec_ctype[:4])

    out_ctype_conv = out_ctype[7] if len(out_ctype) > 4 else out_ctype[:4]
    if CTYPE_TO_PHYSICALTYPE[out_ctype_conv] == 'air wavelength':
        raise NotImplementedError("Conversion to air wavelength is not supported.")
    out_lin_cunit = (LINEAR_CUNIT_DICT[out_ctype_conv] if out_ctype_conv in
                     LINEAR_CUNIT_DICT else outunit)
    out_vcequiv = _parse_velocity_convention(out_ctype_conv)

    # Load the input values
    crval_in = (mywcs.wcs.crval[mywcs.wcs.spec] * inunit)
    # the cdelt matrix may not be correctly populated: need to account for cd,
    # cdelt, and pc
    cdelt_in = (mywcs.pixel_scale_matrix[mywcs.wcs.spec, mywcs.wcs.spec] *
                inunit)

    if in_spec_ctype == 'AWAV':
        warnings.warn("Support for air wavelengths is experimental and only "
                      "works in the forward direction (air->vac, not vac->air).")
        cdelt_in = air_to_vac_deriv(crval_in) * cdelt_in
        crval_in = air_to_vac(crval_in)
        in_spec_ctype = 'WAVE'

    # 1. Convert input to input, linear
    if in_vcequiv is not None and ref_value is not None:
        crval_lin1 = crval_in.to(lin_cunit, u.spectral() + in_vcequiv(ref_value))
    else:
        crval_lin1 = crval_in.to(lin_cunit, u.spectral())
    cdelt_lin1 = cdelt_derivative(crval_in,
                                  cdelt_in,
                                  # equivalent: inunit.physical_type
                                  intype=CTYPE_TO_PHYSICALTYPE[in_spec_ctype[:4]],
                                  outtype=parse_phys_type(lin_cunit),
                                  rest=ref_value,
                                  linear=True)

    # 2. Convert input, linear to output, linear
    if ref_value is None:
        if in_vcequiv is not None:
            pass  # consider raising a ValueError here; not clear if this is valid
        crval_lin2 = crval_lin1.to(out_lin_cunit, u.spectral())
    else:
        # at this stage, the transition can ONLY be relativistic, because the V
        # frame (as a linear frame) is only defined as "apparent velocity"
        crval_lin2 = crval_lin1.to(out_lin_cunit, u.spectral() +
                                   u.doppler_relativistic(ref_value))

    # For cases like VRAD <-> FREQ and VOPT <-> WAVE, this will be linear too:
    linear_middle = in_vcequiv == out_vcequiv

    cdelt_lin2 = cdelt_derivative(crval_lin1, cdelt_lin1,
                                  intype=parse_phys_type(lin_cunit),
                                  outtype=CTYPE_TO_PHYSICALTYPE[out_ctype_conv],
                                  rest=ref_value,
                                  linear=linear_middle)

    # 3. Convert output, linear to output
    if out_vcequiv is not None and ref_value is not None:
        crval_out = crval_lin2.to(outunit, out_vcequiv(ref_value) + u.spectral())
        cdelt_out = cdelt_derivative(crval_lin2,
                                     cdelt_lin2,
                                     intype=CTYPE_TO_PHYSICALTYPE[out_ctype_conv],
                                     outtype=parse_phys_type(outunit),
                                     rest=ref_value,
                                     linear=True).to(outunit)
    else:
        crval_out = crval_lin2.to(outunit, u.spectral())
        cdelt_out = cdelt_lin2.to(outunit, u.spectral())

    if crval_out.unit != cdelt_out.unit:
        # this should not be possible, but it's a sanity check
        raise ValueError("Conversion failed: the units of cdelt and crval don't match.")

    # A cdelt of 0 would be meaningless
    if cdelt_out.value == 0:
        raise ValueError("Conversion failed: the output CDELT would be 0.")

    newwcs = mywcs.deepcopy()
    if hasattr(newwcs.wcs,'cd'):
        newwcs.wcs.cd[newwcs.wcs.spec, newwcs.wcs.spec] = cdelt_out.value
        # todo: would be nice to have an assertion here that no off-diagonal
        # values for the spectral WCS are nonzero, but this is a nontrivial
        # check
    else:
        newwcs.wcs.cdelt[newwcs.wcs.spec] = cdelt_out.value
    newwcs.wcs.cunit[newwcs.wcs.spec] = cdelt_out.unit.to_string(format='fits')
    newwcs.wcs.crval[newwcs.wcs.spec] = crval_out.value
    newwcs.wcs.ctype[newwcs.wcs.spec] = out_ctype
    if rest_value is not None:
        if parse_phys_type(rest_value.unit) == 'frequency':
            newwcs.wcs.restfrq = rest_value.to(u.Hz).value
        elif parse_phys_type(rest_value.unit) == 'length':
            newwcs.wcs.restwav = rest_value.to(u.m).value
        else:
            raise ValueError("Rest Value was specified, but not in frequency or length units")

    newwcs.wcs.set()
    return newwcs


def new_spectral_wcs(wcs, unit, velocity_convention=None, rest_value=None):
    """Returns a new WCS with a different Spectral Axis unit.

    Parameters
    ----------
    wcs : :class:`astropy.wcs.WCS`
        WCS to be converted.
    unit : :class:`~astropy.units.Unit` or str
        Any valid spectral unit: velocity, (wave)length, or frequency.
        Only vacuum units are supported.
    velocity_convention : {'relativistic', 'radio', 'optical', `None`}
        The velocity convention to use for the output velocity axis.
        Required if the output type is velocity. This can be either one
        of the above strings, or an `astropy.units` equivalency.
    rest_value : :class:`~astropy.units.Quantity` or `None`
        A rest wavelength or frequency with appropriate units.  Required if
        output type is velocity. The cube's WCS should include this
        already if the *input* type is velocity, but the WCS's rest
        wavelength/frequency can be overridden with this parameter.

        .. note: This must be the rest frequency/wavelength *in vacuum*,
                 even if your cube has air wavelength units

    Returns
    -------
    newwcs : :class:`astropy.wcs.WCS`
        New WCS in desired unit.

    meta : dict
        Metadata associated with the conversion.

    """
    # Allow string specification of units, for example
    if not isinstance(unit, u.Unit):
        unit = u.Unit(unit)

    # Velocity conventions: required for frq <-> velo
    # convert_spectral_axis will handle the case of no velocity
    # convention specified & one is required
    if velocity_convention in DOPPLER_CONVENTIONS:
        velocity_convention = DOPPLER_CONVENTIONS[velocity_convention]
    elif (velocity_convention is not None and
          velocity_convention not in DOPPLER_CONVENTIONS.values()):
        raise ValueError("Velocity convention must be radio, optical, "
                         "or relativistic.")

    # If rest value is specified, it must be a quantity
    if (rest_value is not None and
            (not hasattr(rest_value, 'unit') or
             not rest_value.unit.is_equivalent(u.m, u.spectral()))):
        raise ValueError("Rest value must be specified as an astropy "
                         "quantity with spectral equivalence.")

    meta = {'Original Unit': wcs.wcs.cunit[wcs.wcs.spec],
            'Original Type': wcs.wcs.ctype[wcs.wcs.spec]}
    out_ctype = ctype_from_vconv(wcs.wcs.ctype[wcs.wcs.spec], unit,
                                 velocity_convention=velocity_convention)
    newwcs = convert_spectral_axis(wcs, unit, out_ctype, rest_value=rest_value)

    return newwcs, meta
