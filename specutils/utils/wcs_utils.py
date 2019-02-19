"""
Utilities for parsing, converting, and manipulating astropy FITS-WCS objects
(i.e., wcsutils objects compliant with Griesen 2006 FITS Paper III)
"""

from astropy import units as u
from astropy import constants
import warnings

def _parse_velocity_convention(vc):
    if vc in (u.doppler_radio, 'radio', 'RADIO', 'VRAD', 'F', 'FREQ'):
        return u.doppler_radio
    elif vc in (u.doppler_optical, 'optical', 'OPTICAL', 'VOPT', 'W', 'WAVE'):
        return u.doppler_optical
    elif vc in (u.doppler_relativistic, 'relativistic', 'RELATIVE', 'VREL',
                'speed', 'V', 'VELO'):
        return u.doppler_relativistic

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
                         'VOPT': 'speed',
                         }

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

def unit_from_header(header):
    """ Retrieve the spectral unit from a header """
    if 'CUNIT3' in header:
        return u.Unit(header['CUNIT3'])

def wcs_unit_scale(unit):
    """
    Determine the appropriate scaling factor to get to the equivalent WCS unit
    """
    for wu in WCS_UNIT_DICT.values():
        if wu.is_equivalent(unit):
            return wu.to(unit)

def determine_vconv_from_ctype(ctype):
    """
    Given a CTYPE, say what velocity convention it is associated with,
    i.e. what unit the velocity is linearly proportional to

    Parameters
    ----------
    ctype : str
        The spectral CTYPE
    """
    if len(ctype) < 5:
        return _parse_velocity_convention(ctype)
    elif len(ctype) == 8:
        return _parse_velocity_convention(ctype[7])
    else:
        raise ValueError("A valid ctype must either have 4 or 8 characters.")

def determine_ctype_from_vconv(ctype, unit, velocity_convention=None):
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
        in_physchar = PHYSICAL_TYPE_TO_CHAR[lin_cunit.physical_type]

    if unit.physical_type == 'speed':
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
        if in_phystype == unit.physical_type:
            # Linear case
            return ALL_CTYPES[in_phystype]
        else:
            # Nonlinear case
            out_physchar = PHYSICAL_TYPE_TO_CTYPE[unit.physical_type]
            return "{type}-{s1}2{s2}".format(type=ALL_CTYPES[unit.physical_type],
                                             s1=in_physchar,
                                             s2=out_physchar)



def get_rest_value_from_wcs(mywcs):
    if mywcs.wcs.restfrq:
        ref_value = mywcs.wcs.restfrq*u.Hz
        return ref_value
    elif mywcs.wcs.restwav:
        ref_value = mywcs.wcs.restwav*u.m
        return ref_value

def convert_spectral_axis(mywcs, outunit, out_ctype, rest_value=None):
    """
    Convert a spectral axis from its unit to a specified out unit with a given output
    ctype

    Only VACUUM units are supported (not air)

    Process:
        1. Convert the input unit to its equivalent linear unit
        2. Convert the input linear unit to the output linear unit
        3. Convert the output linear unit to the output unit
    """

    # If the WCS includes a rest frequency/wavelength, convert it to frequency
    # or wavelength first.  This allows the possibility of changing the rest
    # frequency
    wcs_rv = get_rest_value_from_wcs(mywcs)
    inunit = u.Unit(mywcs.wcs.cunit[mywcs.wcs.spec])
    outunit = u.Unit(outunit)

    # If wcs_rv is set and speed -> speed, then we're changing the reference
    # location and we need to convert to meters or Hz first
    if ((inunit.physical_type == 'speed' and
         outunit.physical_type == 'speed' and
         wcs_rv is not None)):
        mywcs = convert_spectral_axis(mywcs, wcs_rv.unit,
                                      ALL_CTYPES[wcs_rv.unit.physical_type],
                                      rest_value=wcs_rv)
        inunit = u.Unit(mywcs.wcs.cunit[mywcs.wcs.spec])
    elif (inunit.physical_type == 'speed' and outunit.physical_type == 'speed'
          and wcs_rv is None):
        # If there is no reference change, we want an identical WCS, since
        # WCS doesn't know about units *at all*
        newwcs = mywcs.deepcopy()
        return newwcs
        #crval_out = (mywcs.wcs.crval[mywcs.wcs.spec] * inunit).to(outunit)
        #cdelt_out = (mywcs.wcs.cdelt[mywcs.wcs.spec] * inunit).to(outunit)
        #newwcs.wcs.cdelt[newwcs.wcs.spec] = cdelt_out.value
        #newwcs.wcs.cunit[newwcs.wcs.spec] = cdelt_out.unit.to_string(format='fits')
        #newwcs.wcs.crval[newwcs.wcs.spec] = crval_out.value
        #newwcs.wcs.ctype[newwcs.wcs.spec] = out_ctype
        #return newwcs

    in_spec_ctype = mywcs.wcs.ctype[mywcs.wcs.spec]

    # Check whether we need to convert the rest value first
    ref_value = None
    if outunit.physical_type == 'speed':
        if rest_value is None:
            rest_value = wcs_rv
            if rest_value is None:
                raise ValueError("If converting from wavelength/frequency to speed, "
                                 "a reference wavelength/frequency is required.")
        ref_value = rest_value.to(u.Hz, u.spectral())
    elif inunit.physical_type == 'speed':
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
                                  outtype=lin_cunit.physical_type,
                                  rest=ref_value,
                                  linear=True
                                  )

    # 2. Convert input, linear to output, linear
    if ref_value is None:
        if in_vcequiv is not None:
            pass # consider raising a ValueError here; not clear if this is valid
        crval_lin2 = crval_lin1.to(out_lin_cunit, u.spectral())
    else:
        # at this stage, the transition can ONLY be relativistic, because the V
        # frame (as a linear frame) is only defined as "apparent velocity"
        crval_lin2 = crval_lin1.to(out_lin_cunit, u.spectral() +
                                   u.doppler_relativistic(ref_value))

    # For cases like VRAD <-> FREQ and VOPT <-> WAVE, this will be linear too:
    linear_middle = in_vcequiv == out_vcequiv

    cdelt_lin2 = cdelt_derivative(crval_lin1, cdelt_lin1,
                                  intype=lin_cunit.physical_type,
                                  outtype=CTYPE_TO_PHYSICALTYPE[out_ctype_conv],
                                  rest=ref_value,
                                  linear=linear_middle)

    # 3. Convert output, linear to output
    if out_vcequiv is not None and ref_value is not None:
        crval_out = crval_lin2.to(outunit, out_vcequiv(ref_value) + u.spectral())
        #cdelt_out = cdelt_lin2.to(outunit, out_vcequiv(ref_value) + u.spectral())
        cdelt_out = cdelt_derivative(crval_lin2,
                                     cdelt_lin2,
                                     intype=CTYPE_TO_PHYSICALTYPE[out_ctype_conv],
                                     outtype=outunit.physical_type,
                                     rest=ref_value,
                                     linear=True
                                     ).to(outunit)
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
        if rest_value.unit.physical_type == 'frequency':
            newwcs.wcs.restfrq = rest_value.to(u.Hz).value
        elif rest_value.unit.physical_type == 'length':
            newwcs.wcs.restwav = rest_value.to(u.m).value
        else:
            raise ValueError("Rest Value was specified, but not in frequency or length units")

    return newwcs

def cdelt_derivative(crval, cdelt, intype, outtype, linear=False, rest=None):
    if intype == outtype:
        return cdelt
    elif set((outtype,intype)) == set(('length','frequency')):
        # Symmetric equations!
        return (-constants.c / crval**2 * cdelt).to(PHYS_UNIT_DICT[outtype])
    elif outtype in ('frequency','length') and intype == 'speed':
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
    elif outtype == 'speed' and intype in ('frequency','length'):

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


def refraction_index(wavelength, method='Griesen2006', co2=None):
    """
    Calculates the index of refraction of dry air at standard temperature
    and pressure, at different wavelengths, using different methods.

    Parameters
    ----------
    wavelength : `Quantity` object (number or sequence)
        Vacuum wavelengths with an astropy.unit.
    method : str, optional
        Method used to convert wavelengths. Options are:
        'Griesen2006' (default) - from Greisen et al. (2006, A&A 446, 747),
            eqn. 65, standard used by International Unionof Geodesy and Geophysics
        'Edlen1953' - from Edlen (1953, J. Opt. Soc. Am, 43, 339). Standard
            adopted by IAU (resolution No. C15, Commission 44, XXI GA, 1991),
            which refers to Oosterhoff (1957) that uses Edlen (1953). Also used
            by Morton (1991, ApJS, 77, 119), which is frequently cited as IAU source.
        'Edlen1966' - from Edlen (1966, Metrologia 2, 71), rederived constants
            from optical and near UV data.
        'PeckReeder1972' - from Peck & Reeder (1972, J. Opt. Soc. 62), derived
            from additional infrared measurements (up to 1700 nm).
        'Morton2000' - from Morton (2000, ApJS, 130, 403), eqn 8. Used by VALD,
            the Vienna Atomic Line Database. Very similar to Edlen (1966).
        'Ciddor1996' - from Ciddor (1996, Appl. Opt. 35, 1566). Based on
            Peck & Reeder (1972), but updated to account for the changes in
            the international temperature scale and adjust the results for
            CO2 concentration. Arguably most accurate conversion available.
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


def vac_to_air(wavelength, method='Griesen2006', co2=None):
    """
    Converts vacuum to air wavelengths using different methods.

    Parameters
    ----------
    wavelength : `Quantity` object (number or sequence)
        Vacuum wavelengths with an astropy.unit.
    method : str, optional
        One of the methods in refraction_index().
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


def air_to_vac(wavelength, scheme='inversion', method='Griesen2006', co2=None,
               precision=1e-12, maxiter=30):
    """
    Converts air to vacuum wavelengths using different methods.

    Parameters
    ----------
    wavelength : `Quantity` object (number or sequence)
        Air wavelengths with an astropy.unit.
    scheme : str, optional
        How to convert from vacuum to air wavelengths. Options are:
        'inversion' (default) - result is simply the inversion (1 / n) of the
            refraction index of air. Griesen et al. (2006) report that the error
            in naively inverting is less than 10^-9.
        'Piskunov' - uses an analytical solution derived by Nikolai Piskunov
            and used by the Vienna Atomic Line Database (VALD).
        'iteration' - uses an iterative scheme to invert the index of refraction.
    method : str, optional
        Only used if scheme is 'inversion' or 'iteration'. One of the methods
        in refraction_index().
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
                #return wavelength * conv
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
        'Griesen2006' (default) - from Greisen et al. (2006, A&A 446, 747),
            eqn. 66.

    Returns
    -------
    wave_deriv : `Quantity` object (number or sequence)
        Derivative d(wave_vacuum) / d(wave_air).
    """
    assert method.lower() == 'griesen2006', "Only supported method is 'Griesen2006'"
    wlum = wavelength.to(u.um).value
    return (1 + 1e-6 * (287.6155 - 1.62887 / wlum**2 - 0.04080 / wlum**4))
