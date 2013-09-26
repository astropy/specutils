# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Extinction law functions."""

from __future__ import division
import numpy as np
import warnings

def reddening_ccm(wave, ebv=None, a_v=None, r_v=3.1, model='ccm89'):
    """Return the Cardelli, Clayton, Mathis (CCM) reddening values at
    the given wavelength(s).

    Parameters
    ----------
    wave : float or list_like
        Wavelength(s) in Angstroms at which to evaluate the reddening.
    ebv or a_v : float
        E(B-V) differential extinction, or A(V) total V band extinction,
        in magnitudes. Specify exactly one. The two values are related by
        A(V) = R_V * E(B-V).
    r_v : float, optional
        R_V parameter. Default is the standard Milky Way average of 3.1.
    model : {'ccm89', 'gcc09'}, optional
        * 'ccm89' is the default Cardelli, Clayton, & Mathis (1989) [1]_, but
          does include the O'Donnell (1994) parameters to match IDL astrolib.
        * 'gcc09' is Gordon, Cartledge, & Clayton (2009) [2]_. This paper has
          incorrect parameters for the 2175A bump; not yet corrected here.

    Returns
    -------
    reddening : float or `~numpy.ndarray`
        Inverse of flux transmission fraction, equivalent to ``10**(0.4 * a)``
        where ``a`` is the extinction in magnitudes. To deredden spectral flux
        values, multiply by `reddening`. To redden spectral flux values,
        divide by `reddening`.

    Notes
    -----
    Cardelli, Clayton, & Mathis (1989) [1]_ parameterization is used for all
    models. The default parameter values are from CCM except in the optical
    range, where the updated parameters of O'Donnell (1994) [3]_ are used
    (matching the Goddard IDL astrolib routine CCM_UNRED).

    The function is works between 910 A and 3.3 microns, although note the
    default ccm89 model is scientifically valid only at >1250 A.

    Model gcc09 uses the updated UV coefficients of Gordon, Cartledge, & Clayton
    (2009) [2]_, and is valid from 910 A to 3030 A. This function will use CCM89
    at longer wavelengths if GCC09 is selected, but note that the two do not
    connect perfectly smoothly. There is a small discontinuity at 3030 A. Note
    that GCC09 equations 14 and 15 apply to all x>5.9 (the GCC09 paper
    mistakenly states they do not apply at x>8; K. Gordon, priv. comm.).

    References
    ----------
    .. [1] Cardelli, J. A., Clayton, G. C., & Mathis, J. S. 1989, ApJ, 345, 245
    .. [2] Gordon, K. D., Cartledge, S., & Clayton, G. C. 2009, ApJ, 705, 1320
    .. [3] O'Donnell, J. E. 1994, ApJ, 422, 158O

    """

    model = model.lower()
    if model not in ['ccm89','gcc09']:
        raise ValueError('model must be ccm89 or gcc09')
    if (a_v is None) and (ebv is None):
        raise ValueError('Must specify either a_v or ebv')
    if (a_v is not None) and (ebv is not None):
        raise ValueError('Cannot specify both a_v and ebv')
    if a_v is not None:
        ebv = a_v / r_v

    if model == 'gcc09':
        raise ValueError('TEMPORARY: gcc09 currently does 2175A bump '
                         'incorrectly')

    return_scalar = np.isscalar(wave)
    wave = np.asarray(wave)
    x = 1.e4 / np.ravel(wave)  # Inverse microns.

    if np.any(x < 0.3) or np.any(x > 11.):
        raise ValueError('CCM law valid only for wavelengths from '
                         '910 Angstroms to 3.3 microns')
    if np.any(x > 8.) and (model == 'ccm89'):
        warnings.warn('CCM89 should not be used below 1250 A.')
    #    if any(x < 3.3) and any(x > 3.3) and (model == 'gcc09'):
    #        warnings.warn('GCC09 has a discontinuity at 3030 A.')

    a = np.empty_like(x)
    b = np.empty_like(x)

    # NIR
    valid = (0.3 <= x) & (x < 1.1)
    a[valid] = 0.574 * x[valid]**1.61
    b[valid] = -0.527 * x[valid]**1.61

    # optical, using O'Donnell (1994) values
    valid = (1.1 <= x) & (x < 3.3)
    y = x[valid] - 1.82
    coef_a = np.array([-0.505, 1.647, -0.827, -1.718, 1.137, 0.701, -0.609,
                       0.104, 1.])
    coef_b = np.array([3.347, -10.805, 5.491, 11.102, -7.985, -3.989, 2.908,
                       1.952, 0.])
    a[valid] = np.polyval(coef_a, y)
    b[valid] = np.polyval(coef_b, y)

    # UV
    valid = (3.3 <= x) & (x < 8)
    y = x[valid]
    f_a = np.zeros_like(y)
    f_b = np.zeros_like(y)
    select = (y >= 5.9)
    yselect = y[select] - 5.9

    f_a[select] = -0.04473 * yselect**2 - 0.009779 * yselect**3
    f_b[select] = 0.2130 * yselect**2 + 0.1207 * yselect**3
    a[valid] = 1.752 - 0.316*y - (0.104 / ((y-4.67)**2 + 0.341)) + f_a
    b[valid] = -3.090 + 1.825*y + (1.206 / ((y-4.62)**2 + 0.263)) + f_b

    # far-UV CCM89 extrapolation
    valid = (8 <= x) & (x < 11)
    y = x[valid] - 8.
    coef_a = np.array([-0.070, 0.137, -0.628, -1.073])
    coef_b = np.array([0.374, -0.420, 4.257, 13.670])
    a[valid] = np.polyval(coef_a, y)
    b[valid] = np.polyval(coef_b, y)

    # Overwrite UV with GCC09 model if applicable. Not an extrapolation.
    if model == 'gcc09':
        valid = (3.3 <= x) & (x < 11)
        y = x[valid]
        f_a = np.zeros(y.size)
        f_b = np.zeros(y.size)
        select = (5.9 <= y)
        yselect = y[select] - 5.9
        f_a[select] = -0.110 * yselect**2 - 0.0099 * yselect**3
        f_b[select] = 0.537 * yselect**2 + 0.0530 * yselect**3
        a[valid] = 1.896 - 0.372*y - (0.0108 / ((y-4.57)**2 + 0.0422)) + f_a
        b[valid] = -3.503 + 2.057*y + (0.718 / ((y-4.59)**2 + 0.0530*3.1)) + f_b

    a_v = ebv * r_v
    a_lambda = a_v * (a + b / r_v)
    reddening = 10**(0.4 * a_lambda)

    if return_scalar:
        return reddening[0]
    return reddening

def reddening_fm(wave, ebv=None, a_v=None, r_v=3.1, model='f99'):
    """Return the Fitzpatrick & Massa reddening values at
    the given wavelength(s).

    Parameters
    ----------
    wave : float or list_like
        wavelength in Angstroms
    ebv or a_v : float
        E(B-V) differential extinction, or A(V) total V band extinction,
        in magnitudes. Specify exactly one. Related by A(V) = R_V * E(B-V).
    r_v : float, optional
        R_V parameter. Default is the standard Milky Way average of 3.1.
    model : {'f99', 'fm07'}, optional
        * 'f99' is the default Fitzpatrick (1999) [1]_
        * 'fm07' is Fitzpatrick & Massa (2007) [2]_. Currently not R dependent.

    Returns
    -------
    reddening : float or `~numpy.ndarray`
        Inverse of flux transmission fraction, equivalent to ``10**(0.4 * a)``
        where ``a`` is the extinction in magnitudes. To deredden spectral flux
        values, multiply by `reddening`. To redden spectral flux values,
        divide by `reddening`.

    Notes
    -----
    Uses Fitzpatrick (1999) [1]_ by default, which relies on the UV
    parametrization of Fitzpatrick & Massa (1990) [2]_ and spline
    fitting in the optical and IR. This function is defined from 910 A
    to 6 microns, but note the claimed validity goes down only to 1150
    A. The optical spline points are not taken from F99 Table 4, but
    rather updated versions from E. Fitzpatrick (this matches the
    Goddard IDL astrolib routine FM_UNRED).

    The fm07 model uses the Fitzpatrick & Massa (2007) [3]_
    parametrization, which has a slightly different functional
    form. That paper claims it preferable, although it is unclear if
    signficantly (Gordon et al. 2009) [4]_. It is not the literature
    standard, so not default here.

    References
    ----------
    .. [1] Fitzpatrick, E. L. 1999, PASP, 111, 63
    .. [2] Fitpatrick, E. L. & Massa, D. 1990, ApJS, 72, 163
    .. [3] Fitpatrick, E. L. & Massa, D. 2007, ApJ, 663, 320
    .. [4] Gordon, K. D., Cartledge, S., & Clayton, G. C. 2009, ApJ, 705, 1320

    """

    from scipy.interpolate import interp1d

    model = model.lower()
    if model not in ['f99','fm07']:
        raise ValueError('model must be f99 or fm07')
    if (a_v is None) and (ebv is None):
        raise ValueError('Must specify either a_v or ebv')
    if (a_v is not None) and (ebv is not None):
        raise ValueError('Cannot specify both a_v and ebv')
    if a_v is not None:
        ebv = a_v / r_v

    if model == 'fm07':
        raise ValueError('TEMPORARY: fm07 currently not properly R dependent')

    return_scalar = np.isscalar(wave)
    wave = np.asarray(wave)
    x = 1.e4 / np.ravel(wave)  # Inverse microns.

    k = np.zeros_like(x)

    if np.any(x < 0.167) or np.any(x > 11.):
        raise ValueError('fm_dered valid only for wavelengths from 910 A to '
                         '6 microns')

    # UV region
    uvsplit = 1.e4 / 2700.  # Turn 2700A split into inverse microns.
    uv_region = (x >= uvsplit)
    y = x[uv_region]
    k_uv = np.zeros_like(y)

    # Fitzpatrick (1999) model
    if model == 'f99':
        x0, gamma = 4.596, 0.99
        c3, c4 = 3.23, 0.41
        c2 = -0.824 + 4.717 / r_v
        c1 = 2.030 - 3.007 * c2
        D = y**2 / ((y**2 - x0**2)**2 + y**2 * gamma**2)
        F = np.zeros(y.size)
        valid = (y >= 5.9)
        F[valid] = 0.5392 * (y[valid] - 5.9)**2 + 0.05644 * (y[valid] - 5.9)**3
        k_uv = c1 + c2*y + c3*D + c4*F

    # Fitzpatrick & Massa (2007) model
    if model == 'fm07':
        x0, gamma = 4.592, 0.922
        c1, c2, c3, c4, c5 = -0.175, 0.807, 2.991, 0.319, 6.097
        D = y**2 / ((y**2-x0**2)**2 + y**2 * gamma**2)
        valid = (y <= c5)
        k_uv[valid] = c1 + c2*y[valid] + c3*D[valid]
        valid = (y > c5)
        k_uv[valid] = c1 + c2*y[valid] + c3*D[valid] + c4*(y[valid] - c5)**2

    k[uv_region] = k_uv

    # Calculate values for UV spline points to anchor OIR fit
    x_uv_spline = 10000. / np.array([2700., 2600.])
    D = x_uv_spline**2 / ((x_uv_spline**2-x0**2)**2 + x_uv_spline**2 * gamma**2)
    k_uv_spline = c1 + c2*x_uv_spline +c3*D

    # Optical / IR
    OIR_region = (x < uvsplit)
    y = x[OIR_region]
    k_OIR = np.zeros(y.size)

    # Fitzpatrick (1999) model
    if model == 'f99':
        # The OIR anchors are up from IDL astrolib, not F99.
        anchors_extinction = np.array([0, 0.26469*r_v/3.1, 0.82925*r_v/3.1, # IR
            -0.422809 + 1.00270*r_v + 2.13572e-04*r_v**2, # optical
            -5.13540e-02 + 1.00216*r_v - 7.35778e-05*r_v**2,
            0.700127 + 1.00184*r_v - 3.32598e-05*r_v**2,
            (1.19456 + 1.01707*r_v - 5.46959e-03*r_v**2 + 7.97809e-04*r_v**3 +
                -4.45636e-05*r_v**4)])
        anchors_k = np.append(anchors_extinction-r_v, k_uv_spline)
        # Note that interp1d requires that the input abscissa is monotonically
        # _increasing_. This is opposite the usual ordering of a spectrum, but
        # fortunately the _output_ abscissa does not have the same requirement.
        anchors_x = 1e4 / np.array([26500., 12200., 6000., 5470., 4670., 4110.])
        anchors_x = np.append(0., anchors_x)  # For well-behaved spline.
        anchors_x = np.append(anchors_x, x_uv_spline)
        OIR_spline = interp1d(anchors_x, anchors_k, kind='cubic')
        k_OIR = OIR_spline(y)
    # Fitzpatrick & Massa (2007) model
    if model == 'fm07':
        anchors_k_opt = np.array([0., 1.322, 2.055])
        IR_wave = np.array([float('inf'), 4., 2., 1.333, 1.])
        anchors_k_IR = (-0.83 + 0.63*r_v) * IR_wave**-1.84 - r_v
        anchors_k = np.append(anchors_k_IR, anchors_k_opt)
        anchors_k = np.append(anchors_k, k_uv_spline)
        anchors_x = np.array([0., 0.25, 0.50, 0.75, 1.])  # IR
        opt_x = 1.e4 / np.array([5530., 4000., 3300.])  # optical
        anchors_x = np.append(anchors_x, opt_x)
        anchors_x = np.append(anchors_x, x_uv_spline)
        OIR_spline = interp1d(anchors_x, anchors_k, kind='cubic')
        k_OIR = OIR_spline(y)

    k[OIR_region] = k_OIR

    reddening = 10**(0.4 * ebv * (k+r_v))

    if return_scalar:
        return reddening
    return reddening

