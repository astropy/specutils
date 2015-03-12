# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Extinction models.

Classes are callables representing the corresponding extinction
function with a fixed R_V. When calling an extinction function multiple
times with the same R_V, it will be faster to create a class instance with
fixed R_V and use that object to evaluate the extinction. See class
documentation for details.
"""

from __future__ import division
from os import path
import numpy as np

from astropy.io import ascii
from astropy.utils import data as apydata
from astropy import units as u


from specutils import cextinction

try:
    import scipy
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True

__all__ = ['extinction_ccm89', 'extinction_od94', 'extinction_gcc09',
           'extinction_f99', 'extinction_fm07', 'extinction_wd01',
           'extinction_d03', 'extinction', 'reddening',
           'ExtinctionF99', 'ExtinctionD03', 'ExtinctionWD01']

def _process_wave(wave):
    return wave.to(u.angstrom).flatten()

def _process_inputs(wave, ebv, a_v, r_v):
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

def _check_wave(wave, minwave, maxwave):
    if np.any((wave < minwave) | (wave > maxwave)):
        raise ValueError('Wavelengths must be between {0:.2f} and {1:.2f} '
                         'angstroms'.format(minwave, maxwave))
    


def extinction_ccm89(wave, a_v, r_v=3.1):
    """Cardelli, Clayton, & Mathis (1989) extinction model.

    The parameters given in the original paper [1]_ are used.
    The function works between 910 A and 3.3 microns, although note the
    claimed validity is only for wavelength above 1250 A.

    {0}

    Notes
    -----
    In Cardelli, Clayton & Mathis (1989) the mean
    R_V-dependent extinction law, is parameterized as

    .. math::

       <A(\lambda)/A_V> = a(x) + b(x) / R_V

    where the coefficients a(x) and b(x) are functions of
    wavelength. At a wavelength of approximately 5494.5 angstroms (a
    characteristic wavelength for the V band), a(x) = 1 and b(x) = 0,
    so that A(5494.5 angstroms) = A_V. This function returns

    .. math::

       A(\lambda) = A_V (a(x) + b(x) / R_V)

    where A_V can either be specified directly or via E(B-V)
    (by defintion, A_V = R_V * E(B-V)).

    References
    ----------
    .. [1] Cardelli, J. A., Clayton, G. C., & Mathis, J. S. 1989, ApJ, 345, 245
    """

    _check_wave(wave, 909.091 * u.angstrom, 33333.333 * u.angstrom)

    res = cextinction.ccm89(_process_wave(wave).value, a_v, r_v)

    return res.reshape(wave.shape)

def extinction_od94(wave, a_v, r_v=3.1):
    """O'Donnell (1994) extinction model.
    
    Like Cardelli, Clayton, & Mathis (1989) [1]_ but using the O'Donnell
    (1994) [2]_ optical coefficients between 3030 A and 9091 A.

    {0}

    Notes
    -----
    This function matches the Goddard IDL astrolib routine CCM_UNRED.
    From the documentation for that routine:

    1. The CCM curve shows good agreement with the Savage & Mathis (1979)
       [3]_ ultraviolet curve shortward of 1400 A, but is probably
       preferable between 1200 and 1400 A.
    2. Curve is extrapolated between 912 and 1000 A as suggested by
       Longo et al. (1989) [4]_
    3. Valencic et al. (2004) [5]_ revise the ultraviolet CCM
       curve (3.3 -- 8.0 um^-1).    But since their revised curve does
       not connect smoothly with longer and shorter wavelengths, it is
       not included here.

    References
    ----------
    .. [1] Cardelli, J. A., Clayton, G. C., & Mathis, J. S. 1989, ApJ, 345, 245
    .. [2] O'Donnell, J. E. 1994, ApJ, 422, 158O 
    .. [3] Savage & Mathis 1979, ARA&A, 17, 73
    .. [4] Longo et al. 1989, ApJ, 339,474
    .. [5] Valencic et al. 2004, ApJ, 616, 912
    """


    _check_wave(wave, 909.091 * u.angstrom, 33333.333 * u.angstrom)

    res = cextinction.od94(_process_wave(wave).value, a_v, r_v)

    return res.reshape(wave.shape)

def extinction_gcc09(wave, a_v, r_v=3.1):
    """Gordon, Cartledge, & Clayton (2009) extinction model.

    Uses the UV coefficients of Gordon, Cartledge, & Clayton (2009)
    [1]_ between 910 A and 3030 A, otherwise the same as the
    `extinction_od94` function.  Also note that the two do not connect
    perfectly: there is a discontinuity at 3030 A. Note that
    GCC09 equations 14 and 15 apply to all x>5.9 (the GCC09 paper
    mistakenly states they do not apply at x>8; K. Gordon,
    priv. comm.).

    .. warning :: Note that the Gordon, Cartledge, & Clayton (2009) paper
                  has incorrect parameters for the 2175 angstrom bump that
                  have not been corrected here.

    {0}

    References
    ----------
    .. [1] Gordon, K. D., Cartledge, S., & Clayton, G. C. 2009, ApJ, 705, 1320
    """

    _check_wave(wave, 909.091 * u.angstrom, 33333.333 * u.angstrom)

    res = cextinction.gcc09(_process_wave(wave).value, a_v, r_v)

    return res.reshape(wave.shape)

_f99_xknots = 1.e4 / np.array([np.inf, 26500., 12200., 6000., 5470.,
                               4670., 4110., 2700., 2600.])

class ExtinctionF99(object):
    """Fitzpatrick (1999) extinction model with fixed R_V.

    Parameters
    ----------
    r_v : float
        Relation between specific and total extinction, ``a_v = r_v * ebv``.

    Examples
    --------
    Create a callable that gives the extinction law for a given ``r_v``
    and use it:

    >>> f = ExtinctionF99(3.1)
    >>> f(3000., a_v=1.)
    1.7993995521481463

    """

    def __init__(self, a_v, r_v=3.1):
        if not HAS_SCIPY:
            raise ImportError('To use this function scipy needs to be installed')

        from scipy.interpolate import splmake

        self.a_v = a_v
        self.r_v = r_v

        kknots = cextinction.f99kknots(_f99_xknots, self.r_v)

        self._spline = splmake(_f99_xknots, kknots, order=3)

    def __call__(self, wave):
        if not HAS_SCIPY:
            raise ImportError('To use this function scipy needs to be installed')

        from scipy.interpolate import spleval

        wave_shape = wave.shape
        wave = _process_wave(wave)

        _check_wave(wave, 909.091* u.angstrom, 6. * u.micron)

        res = np.empty_like(wave.__array__(), dtype=np.float64)

        # Analytic function in the UV.
        uvmask = wave < (2700. * u.angstrom)
        if np.any(uvmask):
            res[uvmask] = cextinction.f99uv(wave[uvmask].value, self.a_v, self.r_v)

        # Spline in the Optical/IR
        oirmask = ~uvmask
        if np.any(oirmask):
            k = spleval(self._spline, 1. / wave[oirmask].to('micron'))
            res[oirmask] = self.a_v / self.r_v * (k + self.r_v)

        return res.reshape(wave_shape)

def extinction_f99(wave, a_v, r_v=3.1):
    """Fitzpatrick (1999) extinction model.

    Fitzpatrick (1999) [1]_ model which relies on the parametrization
    of Fitzpatrick & Massa (1990) [2]_ in the UV (below 2700 A) and
    spline fitting in the optical and IR. This function is defined
    from 910 A to 6 microns, but note the claimed validity goes down
    only to 1150 A. The optical spline points are not taken from F99
    Table 4, but rather updated versions from E. Fitzpatrick (this
    matches the Goddard IDL astrolib routine FM_UNRED).

    {0}

    References
    ----------
    .. [1] Fitzpatrick, E. L. 1999, PASP, 111, 63
    .. [2] Fitzpatrick, E. L. & Massa, D. 1990, ApJS, 72, 163
    """

    f = ExtinctionF99(a_v, r_v)
    return f(wave)


# fm07 knots for spline
_fm07_r_v = 3.1
_fm07_xknots = np.array([0., 0.25, 0.50, 0.75, 1., 1.e4/5530., 1.e4/4000.,
                        1.e4/3300., 1.e4/2700., 1.e4/2600.])
_fm07_kknots = cextinction.fm07kknots(_fm07_xknots)
try:
    from scipy.interpolate import splmake
    _fm07_spline = splmake(_fm07_xknots, _fm07_kknots, order=3)
except ImportError:
    pass


def extinction_fm07(wave, a_v):
    """Fitzpatrick & Massa (2007) extinction model for R_V = 3.1.

    The Fitzpatrick & Massa (2007) [1]_ model, which has a slightly
    different functional form from that of Fitzpatrick (1999) [3]_
    (`extinction_f99`). Fitzpatrick & Massa (2007) claim it is
    preferable, although it is unclear if signficantly so (Gordon et
    al. 2009 [2]_). Defined from 910 A to 6 microns.

    .. note :: This model is not R_V dependent.

    {0}

    References
    ----------
    .. [1] Fitzpatrick, E. L. & Massa, D. 2007, ApJ, 663, 320
    .. [2] Gordon, K. D., Cartledge, S., & Clayton, G. C. 2009, ApJ, 705, 1320
    .. [3] Fitzpatrick, E. L. 1999, PASP, 111, 63
    """
    if not HAS_SCIPY:
        raise ImportError('To use this function scipy needs to be installed')
    from scipy.interpolate import spleval

    wave_shape = wave.shape
    wave = _process_wave(wave)


    _check_wave(wave, 909.091 * u.angstrom, 6.0 * u.micron)
    res = np.empty_like(wave.__array__(), dtype=np.float64)

    # Simple analytic function in the UV
    uvmask = wave < (2700. * u.angstrom)
    if np.any(uvmask):
        res[uvmask] = cextinction.fm07uv(wave[uvmask].value, a_v)
    
    # Spline in the Optical/IR
    oirmask = ~uvmask
    if np.any(oirmask):
        k = spleval(_fm07_spline, (1. / wave[oirmask].to('micron')).value)
        res[oirmask] = a_v / _fm07_r_v * (k + _fm07_r_v)

    return res.reshape(wave_shape)


prefix = path.join('data', 'extinction_models', 'kext_albedo_WD_MW')
_wd01_fnames = {'3.1': prefix + '_3.1B_60.txt',
                '4.0': prefix + '_4.0B_40.txt',
                '5.5': prefix + '_5.5B_30.txt'}
_d03_fnames = {'3.1': prefix + '_3.1A_60_D03_all.txt',
               '4.0': prefix + '_4.0A_40_D03_all.txt',
               '5.5': prefix + '_5.5A_30_D03_all.txt'}
del prefix

class ExtinctionWD01(object):
    """Weingartner and Draine (2001) extinction model with fixed R_V.

    Parameters
    ----------
    r_v : float
        Relation between specific and total extinction, ``a_v = r_v * ebv``.

    Examples
    --------
    Create a callable that gives the extinction law for a given ``r_v``
    and use it:

    >>> f = ExtinctionWD01(3.1)
    >>> f(3000., a_v=1.)

    Arrays are also accepted and ``ebv`` can be specified instead of ``a_v``:

    >>> f([3000., 4000.], ebv=1./3.1)

    """

    def __init__(self, a_v, r_v):

        if not HAS_SCIPY:
            raise ImportError('To use this function scipy needs to be installed')

        from scipy.interpolate import interp1d

        self.a_v = a_v
        self.r_v = r_v

        fname_key = [item for item in _wd01_fnames.keys() if np.isclose(
            float(item), self.r_v)]

        if len(fname_key) == 0:
            raise ValueError("model only defined for r_v in [3.1, 4.0, 5.5]")
        elif len(fname_key) == 1:
            fname = _wd01_fnames[fname_key[0]]
        else:
            raise ValueError('The given float {0} matches multiple available'
                             ' r_vs [3.1, 4.0, 5.5] - unexpected code error')


        fname = apydata.get_pkg_data_filename(fname)
        data = ascii.read(fname, Reader=ascii.FixedWidth, data_start=51,
                          names=['wave', 'albedo', 'avg_cos', 'C_ext',
                                 'K_abs'],
                          col_starts=[0, 10, 18, 25, 35],
                          col_ends=[9, 17, 24, 34, 42], guess=False)

        # Reverse entries so that they ascend in x (needed for the spline).
        waveknots = np.asarray(data['wave'])[::-1]
        cknots = np.asarray(data['C_ext'])[::-1]
        xknots = 1. / waveknots  # Values in inverse microns.

        # Create a spline just to get normalization.
        spline = interp1d(xknots, cknots)
        cknots = cknots / spline(1.e4 / 5495.)  # Normalize cknots.
        self._spline = interp1d(xknots, cknots)

    def __call__(self, wave):

        wave_shape = wave.shape
        wave = _process_wave(wave)

        x = (1 / wave).to('1/micron')

        res = self.a_v * self._spline(x.value)

        return res.reshape(wave_shape)


def extinction_wd01(wave, a_v, r_v=3.1):
    """Weingartner and Draine (2001) extinction model.

    The Weingartner & Draine (2001) [1]_ dust model.  This model is a
    calculation of the interstellar extinction using a dust model of
    carbonaceous grains and amorphous silicate grains. The
    carbonaceous grains are like PAHs when small and like graphite
    when large. This model is evaluated at discrete wavelengths and
    interpolated between these wavelengths. Grid goes from 1 A to 1000
    microns. The model has been calculated for three different grain
    size distributions which produce interstellar exinctions that look
    like 'ccm89' at Rv = 3.1, Rv = 4.0 and Rv = 5.5.  No interpolation
    to other Rv values is performed, so this model can be evaluated
    only for these values.

    The dust model gives the extinction per H nucleon.  For
    consistency with other extinction laws we normalize this
    extinction law so that it is equal to 1.0 at 5495 angstroms.

    .. note :: Model is not an analytic  function of R_V. Only ``r_v``
               values of 3.1, 4.0 and 5.5 are accepted.

    .. note :: This function reads a table from a file on disk on each call.
               For repeated calls with the same ``r_v``, it will be far faster
               to use the class-based interface `ExtinctionWD01`.

    {0}

    See Also
    --------
    ExtinctionWD01

    References
    ----------
    .. [1] Weingartner, J.C. & Draine, B.T. 2001, ApJ, 548, 296

    """

    f = ExtinctionWD01(a_v, r_v)
    return f(wave)


class ExtinctionD03(ExtinctionWD01):
    """Draine (2003) extinction model with fixed R_V.

    Parameters
    ----------
    r_v : float
        Relation between specific and total extinction, ``a_v = r_v * ebv``.

    Examples
    --------
    Create a callable that gives the extinction law for a given ``r_v``
    and use it:

    >>> f = ExtinctionWD01(3.1)
    >>> f(3000., a_v=1.)

    Arrays are also accepted and ``ebv`` can be specified instead of ``a_v``:

    >>> f([3000., 4000.], ebv=1./3.1)

    """

    def __init__(self, a_v, r_v):
        if not HAS_SCIPY:
            raise ImportError('To use this function scipy needs to be installed')

        from scipy.interpolate import interp1d

        super(ExtinctionD03, self).__init__(a_v, r_v)

        fname_key = [item for item in _wd01_fnames.keys() if np.isclose(
            float(item), self.r_v)]

        if len(fname_key) == 0:
            raise ValueError("model only defined for r_v in [3.1, 4.0, 5.5]")
        elif len(fname_key) == 1:
            fname = _d03_fnames[fname_key[0]]
        else:
            raise ValueError('The given float {0} matches multiple available'
                             ' r_vs [3.1, 4.0, 5.5] - unexpected code error')

        fname = apydata.get_pkg_data_filename(fname)

        data = ascii.read(fname, Reader=ascii.FixedWidth, data_start=67,
                          names=['wave', 'albedo', 'avg_cos', 'C_ext',
                                 'K_abs', 'avg_cos_sq', 'comment'],
                          col_starts=[0, 12, 20, 27, 37, 47, 55],
                          col_ends=[11, 19, 26, 36, 46, 54, 80], guess=False)
        xknots = 1. / np.asarray(data['wave'])
        cknots = np.asarray(data['C_ext'])

        # Create a spline just to get normalization.
        spline = interp1d(xknots, cknots)
        cknots = cknots / spline((1. / (5495. * u.angstrom)).to('1/micron').value)  # Normalize cknots.
        self._spline = interp1d(xknots, cknots)


def extinction_d03(wave, a_v, r_v=3.1):
    """Draine (2003) extinction model.

    The Draine (2003) [2]_ update to WD01 [1]_ where the 
    carbon/PAH abundances relative to 'wd01' have been reduced by a 
    factor of 0.93.

    The dust model gives the extinction per H nucleon.  For
    consistency with other extinction laws we normalize this
    extinction law so that it is equal to 1.0 at 5495 angstroms.

    .. note :: Model is not an analytic  function of R_V. Only ``r_v``
               values of 3.1, 4.0 and  5.5 are accepted.

    {0}

    References
    ----------
    .. [1] Weingartner, J.C. & Draine, B.T. 2001, ApJ, 548, 296
    .. [2] Draine, B.T. 2003, ARA&A, 41, 241
    """
    f = ExtinctionD03(a_v, r_v)
    return f(wave)

_extinction_models = {'ccm89': extinction_ccm89,
                      'od94': extinction_od94,
                      'gcc09': extinction_gcc09,
                      'f99': extinction_f99,
                      'fm07': extinction_fm07,
                      'wd01': extinction_wd01,
                      'd03': extinction_d03}


def extinction(wave, a_v, r_v=3.1, model='od94'):
    """Generic interface for all extinction model functions.

    Parameters
    ----------
    wave : float or list_like
        Wavelength(s) in angstroms.
    a_v : float
        Total V band extinction A(V), in magnitudes. A(V) = R_V * E(B-V).
    r_v : float, optional
        R_V parameter. Default is the standard Milky Way average of 3.1.
    model : {'ccm89', 'od94', 'gcc09', 'f99', 'fm07', 'wd01', d03'}, optional
        Use function ``extinction_[model]``. E.g., for 'ccm89', the function
        ``extinction_ccm89`` is used.

    Returns
    -------
    extinction : float or `~numpy.ndarray`
        Extinction in magnitudes at given wavelengths.

    See Also
    --------
    extinction_ccm89
    extinction_od94
    extinction_gcc09
    extinction_f99
    extinction_fm07
    extinction_wd01
    extinction_d03
    reddening

    Notes
    -----

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

    Examples
    --------

    >>> wave = [2000., 2500., 3000.]
    >>> extinction(wave, a_v=1., r_v=3.1, model='f99')
    array([ 2.76225609,  2.27590036,  1.79939955])

    The extinction scales linearly with ``a_v``. This means
    that when calculating extinction for multiple values of ``a_v``, one can compute extinction ahead of time for a given set of
    wavelengths and then scale by ``a_v`` or ``ebv`` later.  For example:

    >>> a_lambda_over_a_v = extinction(wave, a_v=1.)
    >>> a_v = 0.5
    >>> a_lambda = a_v * a_lambda_over_a_v

    """

    model = model.lower()
    if model not in _extinction_models:
        raise ValueError('unknown model: {0}'.format(model))

    if model == 'fm07':
        if not np.isclose(r_v, 3.1):
            raise ValueError('r_v must be 3.1 for fm07 model')
        return _extinction_models[model](wave, a_v=a_v)
    else:
        return _extinction_models[model](wave, a_v=a_v, r_v=r_v)

def reddening(wave, a_v, r_v=3.1, model='od94'):
    """Inverse of flux transmission fraction at given wavelength(s).

    Parameters
    ----------
    wave : float or list_like
        Wavelength(s) in angstroms at which to evaluate the reddening.
    a_v : float
        Total V band extinction, in magnitudes. A(V) = R_V * E(B-V).
    r_v : float, optional
        R_V parameter. Default is the standard Milky Way average of 3.1.
    model : {'ccm89', 'od94', 'gcc09', 'f99', 'fm07'}, optional

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

    return 10**(0.4 * extinction(wave, a_v, r_v=r_v, model=model))


_func_doc = """

    Parameters
    ----------
    wave : float or list_like
        Wavelength(s) in angstroms.
    a_v : float
        A(V) total V band extinction,
        in magnitudes.
    r_v : float, optional
        R_V parameter. Default is the standard Milky Way average of 3.1.

    Returns
    -------
    extinction : `~numpy.ndarray`
        Extinction in magnitudes at given wavelengths.
    """

for func in [extinction_ccm89, extinction_od94, extinction_gcc09,
             extinction_f99, extinction_fm07, extinction_wd01,
             extinction_d03]:
    func.__doc__ = func.__doc__.format(_func_doc)
