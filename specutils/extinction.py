# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Extinction law functions."""

from __future__ import division
from os import path
import numpy as np
import warnings
from astropy.io import ascii
from astropy.utils import data as apydata

from . import _extinction

__all__ = ['extinction_ccm89', 'extinction_od94', 'extinction_gcc09',
           'extinction_f99', 'extinction', 'reddening']

def process_inputs(wave, ebv, a_v, r_v):
    """Helper function for processing inputs to all extinction_* functions

    Parameters
    ----------
    wave : float or list_like
    ebv, a_v : float or None
    r_v : float

    Returns
    -------
    x : `~numpy.ndarray`, at least 1d
        1e4 / wave
    scalar : bool
        Was input wave a scalar?

    """
    if (a_v is None) and (ebv is None):
        raise ValueError('Must specify either a_v or ebv')
    if (a_v is not None) and (ebv is not None):
        raise ValueError('Cannot specify both a_v and ebv')
    if a_v is None:
        a_v = ebv * r_v

    scalar = np.isscalar(wave)
    wave = np.atleast_1d(wave)
    if wave.ndim > 2:
        raise ValueError("wave cannot be more than 1-d")

    return wave, scalar, a_v

def extinction_ccm89(wave, ebv=None, a_v=None, r_v=3.1):
    wave, scalar, a_v = process_inputs(wave, ebv, a_v, r_v)
    if np.any((wave < 909.09091) | (wave > 33333.33333)):
        raise ValueError('ccm89 law valid only for wavelengths from '
                         '909.09091 to 33333.33333 Angstroms.')
    res = _extinction.ccm89(wave, a_v, r_v)
    if scalar:
        return res[0]
    return res

def extinction_od94(wave, ebv=None, a_v=None, r_v=3.1):
    wave, scalar, a_v = process_inputs(wave, ebv, a_v, r_v)
    if np.any((wave < 909.09091) | (wave > 33333.33333)):
        raise ValueError('od94 law valid only for wavelengths from '
                         '909.09091 to 33333.33333 Angstroms.')
    res = _extinction.od94(wave, a_v, r_v)
    if scalar:
        return res[0]
    return res

def extinction_gcc09(wave, ebv=None, a_v=None, r_v=3.1):
    wave, scalar, a_v = process_inputs(wave, ebv, a_v, r_v)
    if np.any((wave < 909.09091) | (wave > 33333.33333)):
        raise ValueError('gcc09 law valid only for wavelengths from '
                         '909.09091 to 33333.33333 Angstroms.')
    res = _extinction.gcc09(wave, a_v, r_v)
    if scalar:
        return res[0]
    return res

def extinction_f99(wave, ebv=None, a_v=None, r_v=3.1):
    wave, scalar, a_v = process_inputs(wave, ebv, a_v, r_v)
    if np.any((wave < 909.09091) | (wave > 33333.33333)):
        raise ValueError('gcc09 law valid only for wavelengths from '
                         '909.09091 to 33333.33333 Angstroms.')
    res = (a_v / r_v) * (_extinction.f99k(wave, r_v) + r_v)
    if scalar:
        return res[0]
    return res

def _fm07_k_uv(x):
    x0, gamma = 4.592, 0.922
    c1, c2, c3, c4, c5 = -0.175, 0.807, 2.991, 0.319, 6.097
    D = x**2 / ((x**2-x0**2)**2 + x**2 * gamma**2)
    k_uv = np.zeros_like(x)
    valid = (x <= c5)
    k_uv[valid] = c1 + c2*x[valid] + c3*D[valid]
    valid = (x > c5)
    k_uv[valid] = c1 + c2*x[valid] + c3*D[valid] + c4*(x[valid] - c5)**2
    return k_uv

def _f99_like(x, ebv, r_v, model='f99'):
    from scipy.interpolate import interp1d

    if np.any(x < 0.167) or np.any(x > 11.):
        raise ValueError('Wavelength(s) must be between 910 A and 6 um')
    if model == 'fm07' and abs(r_v - 3.1) > 0.001:
        raise ValueError('fm07 model not implementend for r_v != 3.1')

    uvmask = x >= 1.e4 / 2700.
    oirmask = ~uvmask
    k = np.zeros_like(x)

    # UV region
    if model == 'f99':
        k[uvmask] = _f99_k_uv(x[uvmask], r_v)
    if model == 'fm07':
        k[uvmask] = _fm07_k_uv(x[uvmask])

    # Calculate values for UV spline points to anchor OIR fit
    x_uv_spline = 1.e4 / np.array([2700., 2600.])
    d = (x_uv_spline**2 /
         ((x_uv_spline**2 - x0**2)**2 + x_uv_spline**2 * gamma**2))
    k_uv_spline = c1 + c2 * x_uv_spline + c3 * d

    # Optical / IR region
    if model == 'f99':
        anchors_x = 1.e4 / np.array([np.inf, 26500., 12200., 6000., 5470.,
                                     4670., 4110.])

        # The OIR anchors are from IDL astrolib, not F99.
        anchors_extinction = np.array(
            [0.,
             0.26469 * r_v/3.1,  # IR
             0.82925 * r_v/3.1,  # IR
             -0.422809 + 1.00270*r_v + 2.13572e-04*r_v**2,  # optical
             -5.13540e-02 + 1.00216*r_v - 7.35778e-05*r_v**2,
             0.700127 + 1.00184*r_v - 3.32598e-05*r_v**2,
             (1.19456 + 1.01707*r_v - 5.46959e-03*r_v**2 +
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
    k[oirmask] = oir_spline(x[oirmask])

    return ebv * (k + r_v)


d = path.join('data', 'extinction_models')
_wd01_like_fnames = {
    ('wd01', 3.1): path.join(d, 'kext_albedo_WD_MW_3.1B_60.txt'), 
    ('wd01', 4.0): path.join(d, 'kext_albedo_WD_MW_4.0B_40.txt'),
    ('wd01', 5.5): path.join(d, 'kext_albedo_WD_MW_5.5B_30.txt'),
    ('d03', 3.1): path.join(d, 'kext_albedo_WD_MW_3.1A_60_D03_all.txt'),
    ('d03', 4.0): path.join(d, 'kext_albedo_WD_MW_4.0A_40_D03_all.txt'),
    ('d03', 5.5): path.join(d, 'kext_albedo_WD_MW_5.5A_30_D03_all.txt')}
del d
_wd01_like_cache = {}

def _wd01_like(x, ebv, r_v, model=None):
    """
      Read in the dust model, interpolate and normalize

      The dust model gives the extinction per H nucleon.
      For consistency with other extinction laws we 
      normalize this extinction law so that it is equal
      to 1.0 at 5495 angstroms
    """

    from scipy.interpolate import interp1d

    if np.any(x < 0.001) or np.any(x > 1.e4):
        raise ValueError('Wavelength(s) must be between 1 A and 1 mm')

    # If we already read the model in and created a spline, use the cache
    if (model, r_v) in _wd01_like_cache:
        return ebv * _wd01_like_cache[(model, r_v)](x)

    # Get the path to the model file
    try:
        relative_path = _wd01_like_fnames[(model, r_v)]
    except KeyError:
        if model not in ['wd01', 'd03']:
            raise ValueError("model must be 'wd01' or 'd03'")
        else:
            raise ValueError("models only defined for r_v in [3.1, 4.0, 5.5]")
    model_file = apydata.get_pkg_data_filename(relative_path)

    # The d03 and w01 data files have different hard-to-parse formats
    if model == 'd03':
        dust_table = ascii.read(model_file,Reader=ascii.FixedWidth,
                                data_start=67,
                                names=['wave','albedo','avg_cos','C_ext',
                                       'K_abs','avg_cos_sq','comment'],
                                col_starts=[0,12,20,27,37,47,55],
                                col_ends=[11,19,26,36,46,54,80],guess=False)
        waves = np.array(dust_table['wave'])
        c_vals = np.array(dust_table['C_ext'])

    if model == 'wd01':
        dust_table = ascii.read(model_file,Reader=ascii.FixedWidth,
                                data_start=51,
                                names=['wave','albedo','avg_cos','C_ext',
                                       'K_abs'],
                                col_starts=[0,10,18,25,35],
                                col_ends=  [9,17,24,34,42],guess=False)

        #We have to reverse the entries for wd01
        waves = np.array(dust_table['wave'])[::-1]
        c_vals = np.array(dust_table['C_ext'])[::-1]

    #Put waves into inverse microns
    waves_ang = 1./waves

    # Create a spline to get normalization at effective V band
    ext_spline = interp1d(waves_ang,c_vals)
    normalization = ext_spline(1.e4/5495.)
    
    # Normalize to ebv=1. and create a new spline with the normalized values
    c_vals /= normalization * r_v
    ext_spline = interp1d(waves_ang, c_vals)
    _wd01_like_cache[(model, r_v)] = ext_spline  # cache the spline
    
    return ebv * ext_spline(x)



_extinction_models = {'ccm89': extinction_ccm89,
                      'od94': extinction_od94,
                      'gcc09': extinction_gcc09,
                      'f99': extinction_f99}

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
    model : {'ccm89', 'od94', 'gcc09', 'f99', 'fm07', 'wd01', d03'}, optional
        * 'ccm89': Cardelli, Clayton, & Mathis (1989).
        * 'od94': Like 'ccm89' but using the O'Donnell (1994) optical
          coefficients between 3030 A and 9091 A.
        * 'gcc09': Gordon, Cartledge, & Clayton (2009) model at
          wavelengths below 3030 A, otherwise the same as O'Donnell (1994).
          Note that the GCC09 paper has incorrect parameters for the
          2175A bump, and these incorrect parameters **are** used here.
        * 'f99': Fitzpatrick (1999) model. **[DEFAULT]**
        * 'fm07': Fitzpatrick & Massa (2007) model. This model is not
          currently not R dependent, so can only be used with ``r_v=3.1``.
        * 'wd01': Weingartner and Draine (2001) dust model.
        * 'd03': Draine (2003) model. Like the 'wd01' but with lower
          PAH C abundance relative to H.

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
      the Goddard IDL astrolib routine FM_UNRED). This model is the
      DEFAULT extinction law used.

    * **'fm07'** The Fitzpatrick & Massa (2007) [6]_ model, which has
      a slightly different functional form from 'f99'. Fitzpatrick &
      Massa (2007) claim it is preferable, although it is unclear if
      signficantly so (Gordon et al. 2009). Defined from 910 A to 6
      microns.

    * **'wd01'** The Weingartner & Draine (2001) [10]_ dust model.
      This model is a calculation of the interstellar extinction
      using a dust model of carbonaceous grains and amorphous 
      silicate grains. The carbonaceous grains are like PAHs when 
      small and like graphite when large. This model is evaluated
      at discrete wavelengths and interpolated between these 
      wavelengths. Grid goes from 1 A to 1000 microns. The model 
      has been calculated for three different grain size 
      distributions which produce interstellar exinctions that 
      look like 'ccm89' at Rv = 3.1, Rv = 4.0 and Rv = 5.5.
      No interpolation to other Rv values is performed, so this 
      model can be evaluated only for these values.

    * **'d03'** The Draine (2003) [11]_ update to 'wd01' where the 
      carbon/PAH abundances relative to 'wd01' have been reduced by a 
      factor of 0.93. 

    .. warning :: The 'gcc09' model currently does the 2175 Angstrom bump
                  incorrectly.

    .. warning :: The 'fm07' model is not properly R dependent. A ValueError
                  is raised if r_v values other than 3.1 are used.

    .. warning :: The 'wd01' and 'd03' models are not analytic 
                  functions of R. Only r_v = 3.1, r_v = 4.0 and
                  r_v = 5.5 are accepted. Otherwise a ValueError
                  is raised.


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

       models = ['ccm89', 'od94', 'gcc09', 'f99', 'fm07','wd01','d03']
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
       plt.ylim(ymin=-0.4,ymax=0.4)
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
    .. [10] Weingartner, J.C. & Draine, B.T. 2001, ApJ, 548, 296
    .. [11] Draine, B.T. 2003, ARA&A, 41, 241

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
    if model not in _extinction_models:
        raise ValueError('unknown model: {0}'.format(model))
    return _extinction_models[model](wave, ebv=ebv, a_v=a_v, r_v=r_v) 

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
        * 'wd01': Weingartner & Draine (2001)
        * 'd03': Draine (2003)

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
