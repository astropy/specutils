# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Extinction law functions."""

from __future__ import division
import numpy as np
import warnings

def extinction(wave, ebv=None, a_v=None, r_v=3.1, model='f99'):
    """Return extinction in magnitudes at given wavelength(s).

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
    model : {'ccm89', 'od94', 'gcc09', 'f99', 'fm07'}, optional
        * 'ccm89': Cardelli, Clayton, & Mathis (1989).
        * 'od94': Like 'ccm89' but using the O'Donnell (1994) optical
          coefficients between 3030 A and 9091 A.
        * 'gcc09': Gordon, Cartledge, & Clayton (2009) model at
          wavelengths below 3030 A, otherwise the same as O'Donnell (1994).
          Note that the GCC09 paper has incorrect parameters for the
          2175A bump, and these incorrect parameters **are** used here.
        * 'f99': Fitzpatrick (1999) model.
        * 'fm07': Fitzpatrick & Massa (2007) model. This model is not
          currently not R dependent, so can only be used with ``r_v=3.1``.

    Returns
    -------
    extinction : float or `~numpy.ndarray`
        Extinction in magnitudes at given wavelengths.

    See Also
    --------
    reddening

    Notes
    -----

    Description of model options:

    * **'ccm89'** The Cardelli, Clayton, & Mathis (1989) [1]_
      parameterization with coefficients given in that paper. The
      function works between 910 A and 3.3 microns, although note the
      model is scientifically valid only at >1250 A. A warning is
      issued if any values are between 910 A and 1250 A.

    * **'od94'** Same as 'ccm89' but uses optical coefficients from
      O'Donnell (1994) [3]_. This matches the Goddard IDL astrolib routine
      CCM_UNRED.

    * **'gcc09'** uses the updated UV coefficients of Gordon,
      Cartledge, & Clayton (2009) [2]_, and is valid from 910 A to
      3030 A. This function will use the 'od94' model at longer
      wavelengths, but note that the two do not connect perfectly
      smoothly: there is a small discontinuity at 3030 A. Note that
      GCC09 equations 14 and 15 apply to all x>5.9 (the GCC09 paper
      mistakenly states they do not apply at x>8; K. Gordon,
      priv. comm.).

    * **'f99'** Fitzpatrick (1999) [4]_ model which relies on the UV
      parametrization of Fitzpatrick & Massa (1990) [5]_ and spline
      fitting in the optical and IR. This function is defined from 910
      A to 6 microns, but note the claimed validity goes down only to
      1150 A. The optical spline points are not taken from F99 Table
      4, but rather updated versions from E. Fitzpatrick (this matches
      the Goddard IDL astrolib routine FM_UNRED).

    * **'fm07'** The Fitzpatrick & Massa (2007) [6]_ model, which has
      a slightly different functional form from 'f99'. Fitzpatrick &
      Massa (2007) claim it is preferable, although it is unclear if
      signficantly so (Gordon et al. 2009). Defined from 910 A to 6
      microns.

    .. warning :: The 'gcc09' model currently does the 2175 Angstrom bump
                  incorrectly.

    .. warning :: The 'fm07' model is not properly R dependent. A ValueError
                  is raised if r_v values other than 3.1 are used.

    **ccm89, od94 & gcc09 models**

    In Cardelli, Clayton & Mathis (1989) the mean
    R_V-dependent extinction law, is parameterized as

    .. math::

       <A(\lambda)/A_V> = a(x) + b(x) / R_V

    where the coefficients a(x) and b(x) are functions of
    wavelength. At a wavelength of approximately 5494.5 Angstroms (a
    characteristic wavelength for the V band), a(x) = 1 and b(x) = 0,
    so that A(5494.5 Angstroms) = A_V. This function returns

    .. math::

       A(\lambda) = A_V (a(x) + b(x) / R_V)

    where A_V can either be specified directly or via E(B-V)
    (by defintion, A_V = R_V * E(B-V)).

    Notes from the IDL routine CCM_UNRED:

    1. The CCM curve shows good agreement with the Savage & Mathis (1979)
       [7]_ ultraviolet curve shortward of 1400 A, but is probably
       preferable between 1200 and 1400 A.
    2. Curve is extrapolated between 912 and 1000 A as suggested by
       Longo et al. (1989) [8]_
    3. Valencic et al. (2004) [9]_ revise the ultraviolet CCM
       curve (3.3 -- 8.0 um^{-1}).    But since their revised curve does
       not connect smoothly with longer and shorter wavelengths, it is
       not included here.

    
    **Visual comparison of models**

    The plot below shows a comparison of the models for
    ``r_v=3.1``. The shaded regions show the limits of claimed
    validity for the f99 model (> 1150 A) and the ccm89/od94 models (>
    1250 A). The vertical dotted lines indicate transition wavelengths in
    the models (3030.3 A and 9090.9 A for ccm89, od94 and gcc09;
    2700. A for f99, fm07).

    .. plot::

       import numpy as np
       import matplotlib.pyplot as plt
       from mpl_toolkits.axes_grid1 import make_axes_locatable
       from specutils.extinction import extinction

       models = ['ccm89', 'od94', 'gcc09', 'f99', 'fm07']
       wave = np.logspace(np.log10(910.), np.log10(30000.), 2000)
       a_lambda = {model: extinction(wave, a_v=1., model=model)
                   for model in models}

       fig = plt.figure(figsize=(8.5, 6.))

       ax = plt.axes()
       for model in models:
           plt.plot(wave, a_lambda[model], label=model)
       plt.axvline(x=2700., ls=':', c='k')
       plt.axvline(x=3030.3030, ls=':', c='k')
       plt.axvline(x=9090.9091, ls=':', c='k')
       plt.axvspan(wave[0], 1150., fc='0.8', ec='none', zorder=-1000)
       plt.axvspan(1150., 1250., fc='0.9', ec='none', zorder=-1000)    
       plt.text(0.67, 0.95, '$R_V = 3.1$', transform=ax.transAxes, va='top',
                size='x-large')
       plt.ylabel('Extinction ($A(\lambda)$ / $A_V$)')
       plt.legend()
       plt.setp(ax.get_xticklabels(), visible=False)

       divider = make_axes_locatable(ax)
       axresid = divider.append_axes("bottom", size=2.0, pad=0.2, sharex=ax)
       for model in models:
           plt.plot(wave, a_lambda[model] - a_lambda['f99'])
       plt.axvline(x=2700., ls=':', c='k')
       plt.axvline(x=3030.3030, ls=':', c='k')
       plt.axvline(x=9090.9091, ls=':', c='k')
       plt.axvspan(wave[0], 1150., fc='0.8', ec='none', zorder=-1000)
       plt.axvspan(1150., 1250., fc='0.9', ec='none', zorder=-1000)
       plt.xlim(wave[0], wave[-1])
       plt.ylim(ymax=0.4)
       plt.ylabel('residual from f99')
       plt.xlabel('Wavelength ($\AA$)')

       ax.set_xscale('log')
       axresid.set_xscale('log')
       plt.tight_layout()

       fig.show()

    References
    ----------
    .. [1] Cardelli, J. A., Clayton, G. C., & Mathis, J. S. 1989, ApJ, 345, 245
    .. [2] Gordon, K. D., Cartledge, S., & Clayton, G. C. 2009, ApJ, 705, 1320
    .. [3] O'Donnell, J. E. 1994, ApJ, 422, 158O
    .. [4] Fitzpatrick, E. L. 1999, PASP, 111, 63
    .. [5] Fitpatrick, E. L. & Massa, D. 1990, ApJS, 72, 163
    .. [6] Fitpatrick, E. L. & Massa, D. 2007, ApJ, 663, 320
    .. [7] Savage & Mathis 1979, ARA&A, 17, 73
    .. [8] Longo et al. 1989, ApJ, 339,474
    .. [9] Valencic et al. 2004, ApJ, 616, 912

    Examples
    --------

    >>> import numpy as np
    >>> from specutils.extinction import extinction
    >>> wave = np.array([2000., 2500., 3000.])
    >>> extinction(wave, a_v=1., r_v=3.1, model='f99')
    array([ 2.76225609,  2.27590036,  1.79939955])

    The extinction scales linearly with `a_v` or `ebv`. This means
    that when calculating extinction for multiple values of `a_v` or
    `ebv`, one can compute extinction ahead of time for a given set of
    wavelengths and then scale by `a_v` or `ebv` later.  For example:

    >>> a_lambda_over_a_v = extinction(wave, a_v=1.)  # somewhat slow
    >>> a_v = 0.5
    >>> a_lambda = a_v * a_lambda_over_a_v  # relatively fast

    Similarly for `ebv`:

    >>> a_lambda_over_ebv = extinction(wave, ebv=1.)  # somewhat slow
    >>> ebv = 0.1
    >>> a_lambda = ebv * a_lambda_over_ebv  # relatively fast

    """

    model = model.lower()
    if (a_v is None) and (ebv is None):
        raise ValueError('Must specify either a_v or ebv')
    if (a_v is not None) and (ebv is not None):
        raise ValueError('Cannot specify both a_v and ebv')
    if a_v is not None:
        ebv = a_v / r_v

    return_scalar = np.isscalar(wave)
    wave = np.asarray(wave)
    x = 1.e4 / np.ravel(wave)  # Inverse microns.

    if model == 'ccm89':
        a_lambda = _ccm89_like(x, ebv, r_v, ccm89_coeffs_a, ccm89_coeffs_b)
    elif model == 'od94':
        a_lambda = _ccm89_like(x, ebv, r_v, od94_coeffs_a, od94_coeffs_b)
    elif model == 'gcc09':
        uv = (x >= 3.3) & (x < 11)
        non_uv = ~uv
        if np.any(uv) and np.any(non_uv):
            warnings.warn('gcc09 model discontinuous at 3030 A.')
        a_lambda = np.empty_like(x)
        a_lambda[uv] = _gcc09(x[uv], ebv, r_v)
        a_lambda[non_uv] = _ccm89_like(x[non_uv], ebv, r_v, od94_coeffs_a,
                                       od94_coeffs_b)
    elif model == 'f99':
        a_lambda = _f99_like(x, ebv, r_v, model='f99')
    elif model == 'fm07':
        a_lambda = _f99_like(x, ebv, r_v, model='fm07')
    else:
        raise ValueError('unknown model: {}'.format(model))

    if return_scalar:
        return a_lambda[0]
    return a_lambda

def reddening(wave, ebv=None, a_v=None, r_v=3.1, model='od94'):
    """Return reddening (inverse of flux transmission fraction) at given
    wavelength(s).

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
    model : {'ccm89', 'od94', 'gcc09', 'f99', 'fm07'}, optional
        * 'ccm89': Cardelli, Clayton, & Mathis (1989).
        * 'od94': Like 'ccm89' but using the O'Donnell (1994) optical
          coefficients.
        * 'gcc09': Gordon, Cartledge, & Clayton (2009) model at
          wavelengths below 3030A, otherwise the same as O'Donnell (1994).
          This paper has incorrect parameters for the 2175A bump; not yet
          corrected here.
        * 'f99': Fitzpatrick (1999).
        * 'fm07': Fitzpatrick & Massa (2007). Currently not R dependent.

    Returns
    -------
    reddening : float or `~numpy.ndarray`
        Inverse of flux transmission fraction, equivalent to
        ``10**(0.4 * extinction(wave))``. To deredden spectra,
        multiply flux values by these value(s). To redden spectra, divide
        flux values by these value(s).

    See Also
    --------
    extinction

    """

    return 10**(0.4 * extinction(wave, ebv=ebv, a_v=a_v, r_v=r_v, model=model))

def _gcc09(x, ebv, r_v):
    f_a = np.zeros_like(x)
    f_b = np.zeros_like(x)
    select = x >= 5.9
    y = x[select] - 5.9
    f_a[select] = -0.110 * y**2 - 0.0099 * y**3
    f_b[select] = 0.537 * y**2 + 0.0530 * y**3
    a = 1.896 - 0.372 * x - (0.0108 / ((x - 4.57)**2 + 0.0422)) + f_a
    b = -3.503 + 2.057 * x + (0.718 / ((x - 4.59)**2 + 0.0530 * 3.1)) + f_b

    return  ebv * (r_v * a + b)

# Optical/NIR coefficientsfor ccm89 and od94 models.
# These are ordered like: c[0] * x^7 + c[1] * x^6 + ... + c[7] * x^0
# corresponding to numpy.polyval() usage.
ccm89_coeffs_a = np.array([0.32999, -0.77530, 0.01979, 0.72085, -0.02427,
                           -0.50447, 0.17699, 1.])
ccm89_coeffs_b = np.array([-2.09002, 5.30260, -0.62251, -5.38434, 1.07233,
                            2.28305, 1.41338, 0.])
od94_coeffs_a = np.array([-0.505, 1.647, -0.827, -1.718, 1.137, 0.701, -0.609,
                           0.104, 1.])
od94_coeffs_b = np.array([3.347, -10.805, 5.491, 11.102, -7.985, -3.989, 2.908,
                          1.952, 0.])

def _ccm89_like(x, ebv, r_v, optical_coeffs_a, optical_coeffs_b):
    if np.any(x < 0.3) or np.any(x > 11.):
        raise ValueError('CCM law valid only for wavelengths from '
                         '910 Angstroms to 3.3 microns')

    a = np.empty_like(x)
    b = np.empty_like(x)

    # Near Infrared.
    valid = (0.3 <= x) & (x < 1.1)
    a[valid] = 0.574 * x[valid]**1.61
    b[valid] = -0.527 * x[valid]**1.61

    # Optical.
    valid = (1.1 <= x) & (x < 3.3)
    y = x[valid] - 1.82
    a[valid] = np.polyval(optical_coeffs_a, y)
    b[valid] = np.polyval(optical_coeffs_b, y)

    # Ultraviolet.
    valid = (3.3 <= x) & (x < 8.)
    y = x[valid]
    f_a = np.zeros_like(y)
    f_b = np.zeros_like(y)
    select = (y >= 5.9)
    yselect = y[select] - 5.9
    f_a[select] = -0.04473 * yselect**2 - 0.009779 * yselect**3
    f_b[select] = 0.2130 * yselect**2 + 0.1207 * yselect**3
    a[valid] = 1.752 - 0.316*y - (0.104 / ((y-4.67)**2 + 0.341)) + f_a
    b[valid] = -3.090 + 1.825*y + (1.206 / ((y-4.62)**2 + 0.263)) + f_b

    # Far-UV (CCM89 extrapolation)
    valid = (8. <= x) & (x < 11.)
    if np.any(valid):
        warnings.warn('ccm89 and od94 models should not be used below 1250 A.')
    y = x[valid] - 8.
    coef_a = np.array([-0.070, 0.137, -0.628, -1.073])
    coef_b = np.array([0.374, -0.420, 4.257, 13.670])
    a[valid] = np.polyval(coef_a, y)
    b[valid] = np.polyval(coef_b, y)

    return ebv * (r_v * a + b)

def _f99_like(x, ebv, r_v, model='f99'):
    from scipy.interpolate import interp1d

    if np.any(x < 0.167) or np.any(x > 11.):
        raise ValueError('Wavelength(s) must be between 910 A and 6 um')
    if model == 'fm07' and abs(r_v - 3.1) > 0.001:
        raise ValueError('fm07 model not implementend for r_v != 3.1')

    k = np.zeros_like(x)
    uv_region = (x >= 1.e4 / 2700.)
    oir_region = ~uv_region

    # UV region
    y = x[uv_region]
    if model == 'f99':
        x0, gamma = 4.596, 0.99
        c3, c4, c5 = 3.23, 0.41, 5.9
        c2 = -0.824 + 4.717 / r_v
        c1 = 2.030 - 3.007 * c2
        d = y**2 / ((y**2 - x0**2)**2 + y**2 * gamma**2)
        f = np.zeros_like(y)
        valid = (y >= c5)
        f[valid] = 0.5392 * (y[valid] - c5)**2 + 0.05644 * (y[valid] - c5)**3
        k_uv = c1 + c2 * y + c3 * d + c4 * f
    if model == 'fm07':
        x0, gamma = 4.592, 0.922
        c1, c2, c3, c4, c5 = -0.175, 0.807, 2.991, 0.319, 6.097
        D = y**2 / ((y**2-x0**2)**2 + y**2 * gamma**2)
        k_uv = np.zeros_like(y)
        valid = (y <= c5)
        k_uv[valid] = c1 + c2*y[valid] + c3*D[valid]
        valid = (y > c5)
        k_uv[valid] = c1 + c2*y[valid] + c3*D[valid] + c4*(y[valid] - c5)**2
    k[uv_region] = k_uv

    # Calculate values for UV spline points to anchor OIR fit
    x_uv_spline = 1.e4 / np.array([2700., 2600.])
    d = (x_uv_spline**2 /
         ((x_uv_spline**2 - x0**2)**2 + x_uv_spline**2 * gamma**2))
    k_uv_spline = c1 + c2 * x_uv_spline + c3 * d

    # Optical / IR region
    y = x[oir_region]
    if model == 'f99':
        anchors_x = 1.e4 / np.array([np.inf, 26500., 12200., 6000., 5470.,
                                     4670., 4110.])

        # The OIR anchors are from IDL astrolib, not F99.
        anchors_extinction = np.array(
            [0.,
             0.26469 * r_v / 3.1,  # IR
             0.82925 * r_v / 3.1,  # IR
             -0.422809 + 1.00270 * r_v + 2.13572e-04 * r_v**2,  # optical
             -5.13540e-02 + 1.00216 * r_v - 7.35778e-05 * r_v**2,
             0.700127 + 1.00184 * r_v - 3.32598e-05 * r_v**2,
             (1.19456 + 1.01707 * r_v - 5.46959e-03 * r_v**2 +
              7.97809e-04 * r_v**3 - 4.45636e-05 * r_v**4)]
            )

        anchors_x = np.append(anchors_x, x_uv_spline)
        anchors_k = np.append(anchors_extinction - r_v, k_uv_spline)

    if model == 'fm07':
        anchors_x_ir = np.array([0., 0.25, 0.50, 0.75, 1.])
        anchors_k_ir = (-0.83 + 0.63 * r_v) * anchors_x_ir**1.84 - r_v
        anchors_x_opt = np.array([5530., 4000., 3300.])
        anchors_k_opt = np.array([0., 1.322, 2.055])

        anchors_x = np.append(anchors_x_ir, anchors_x_opt)
        anchors_k = np.append(anchors_k_ir, anchors_k_opt)

        anchors_x = np.append(anchors_x, x_uv_spline)
        anchors_k = np.append(anchors_k, k_uv_spline)
    
    # Note that interp1d requires that the input abscissa is monotonically
    # _increasing_. This is opposite the usual ordering of a spectrum, but
    # fortunately the _output_ abscissa does not have the same requirement.
    oir_spline = interp1d(anchors_x, anchors_k, kind='cubic')
    k[oir_region] = oir_spline(y)

    return ebv * (k + r_v)
