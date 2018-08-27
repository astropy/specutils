from __future__ import print_function, division, absolute_import

from astropy import modeling
from astropy.modeling import models, fitting
from astropy.nddata import StdDevUncertainty
from ..spectra import Spectrum1D

import numpy as np
from scipy import interpolate

import logging
import warnings

__all__ = ['fit_continuum_generic', 'fit_continuum_linetools']

def fit_continuum_generic(spectrum,
                          model=None, fitter=None,
                          sigma=3.0, sigma_lower=None, sigma_upper=None, iters=5,
                          exclude_regions=[],
                          full_output=False):
    """
    Fit a generic continuum model to a spectrum.
    
    The default algorithm is iterative sigma clipping
    
    Parameters
    ----------
    spectrum : `~specutils.Spectrum1D`
        The `~specutils.Spectrum1D` object to which a continuum model is fit

    model : `XXXX`
        The type of model to use for the continuum.
        astropy.modeling.models
        Must either be astropy.modeling.Fittable1DModel
        or the string "spline" (since this is not currently implemented)
        Default: models.Chebyshev1D(3)
        
    fitter : `XXXX`
        The type of fitter to use for the continuum.
        astropy.modeling.fitting
        Default: fitting.LevMarLSQFitter()
        
    sigma : float, optional
        The number of standard deviations to use for both lower and upper clipping limit.
        Defaults to 3.0
        
    sigma_lower : float or None, optional
        Number of standard deviations for lower bound clipping limit.
        If None (default), then `sigma` is used.
        
    sigma_upper : float or None, optional
        Number of standard deviations for upper bound clipping limit.
        If None (default), then `sigma` is used.
        
    iters : int or None, optional
        Number of iterations to perform sigma clipping.
        If None, clips until convergence achieved.
        Defaults to 5
        
    exclude_regions : list of tuples, optional
        A list of dispersion regions to exclude.
        Each tuple must be sorted.
        e.g. [(6555,6575)]

    full_output : bool, optional
        If True, return more information.
        Currently, just the model and the pixels-used boolean array
        
    Returns
    -------
    continuum_model : `XXXX`
        Output `XXXX` which is a model for the continuum

    Raises
    ------
    ValueError
       In the case that ``spectrum`` .... is not the correct type
    
    """
    
    ## Parameter checks
    if not isinstance(spectrum, Spectrum1D):
        raise ValueError('The spectrum parameter must be a Spectrum1D object')
    for exclude_region in exclude_regions:
        if len(exclude_region) != 2:
            raise ValueError('All exclusion regions must be of length 2')
        if exclude_region[0] >= exclude_region[1]:
            raise ValueError('All exclusion regions must be (low, high)')
    
    ## Set default model and fitter
    if model is None:
        logging.info("Using Chebyshev1D(3) as default continuum model")
        model = models.Chebyshev1D(3)
    if fitter is None:
        fitter = fitting.LevMarLSQFitter()
    if not isinstance(model, modeling.FittableModel):
        raise ValueError('The model parameter must be a astropy.modeling.FittableModel object')
    ## TODO this is waiting on a refactor in fitting to work
    #if not isinstance(fitter, fitting.Fitter):
    #    raise ValueError('The model parameter must be a astropy.modeling.fitting.Fitter object')

    ## Get input spectrum data
    x = spectrum.spectral_axis.value
    y = spectrum.flux.value
    
    ## Set up valid pixels mask
    ## Exclude non-finite values
    good = np.isfinite(y)
    ## Exclude regions
    for (excl1, excl2) in exclude_regions:
        good[np.logical_and(x > excl1, x < excl2)] = False
    
    ## Set up sigma clipping
    if sigma_lower is None: sigma_lower = sigma
    if sigma_upper is None: sigma_upper = sigma

    for i_iter in range(iters):
        logging.info("Iter {}: Fitting {}/{} pixels".format(i_iter, good.sum(), len(good)))
        ## Fit model
        ## TODO include data uncertainties
        new_model = fitter(model, x[good], y[good])
        
        ## Sigma clip
        difference = new_model(x) - y
        finite = np.isfinite(difference)
        sigma_difference = difference / np.std(difference[np.logical_and(good, finite)])
        good[sigma_difference > sigma_upper] = False
        good[sigma_difference < -sigma_lower] = False
        
    model = new_model
    if full_output:
        return model, good
    return model

def fit_continuum_linetools(spec, edges=None, ax=None, debug=False, kind="QSO", **kwargs):
    """
    A direct port of the linetools continuum normalization algorithm by X Prochaska
    https://github.com/linetools/linetools/blob/master/linetools/analysis/continuum.py
    
    The only changes are switching to Scipy's Akima1D interpolator and changing the relevant syntax
    """
    assert kind in ["QSO"], kind
    if not isinstance(spec, Spectrum1D):
        raise ValueError('The spectrum parameter must be a Spectrum1D object')
    
    ### To start, we define all the functions here to avoid namespace bloat, but this can be fixed later
    ### The goal is to have the same algorithm but with flexible wavelength chunks for other object types
    
    def make_chunks_qso(wa, redshift, divmult=1, forest_divmult=1, debug=False):
        """ Generate a series of wavelength chunks for use by
        prepare_knots, assuming a QSO spectrum.
        """
    
        cond = np.isnan(wa)
        if np.any(cond):
            warnings.warn('Some wavelengths are NaN, ignoring these pixels.')
            wa = wa[~cond]
            assert len(wa) > 0
    
        zp1 = 1 + redshift
        div = np.rec.fromrecords([(200. , 500. , 25),
                                  (500. , 800. , 25),
                                  (800. , 1190., 25),
                                  (1190., 1213.,  4),
                                  (1213., 1230.,  6),
                                  (1230., 1263.,  6),
                                  (1263., 1290.,  5),
                                  (1290., 1340.,  5),
                                  (1340., 1370.,  2),
                                  (1370., 1410.,  5),
                                  (1410., 1515.,  5),
                                  (1515., 1600., 15),
                                  (1600., 1800.,  8),
                                  (1800., 1900.,  5),
                                  (1900., 1940.,  5),
                                  (1940., 2240., 15),
                                  (2240., 3000., 25),
                                  (3000., 6000., 80),
                                  (6000., 20000., 100),
                                  ], names=str('left,right,num'))
    
        div.num[2:] = np.ceil(div.num[2:] * divmult)
        div.num[:2] = np.ceil(div.num[:2] * forest_divmult)
        div.left *= zp1
        div.right *= zp1
        if debug:
            print(div.tolist())
        temp = [np.linspace(left, right, n+1)[:-1] for left,right,n in div]
        edges = np.concatenate(temp)
    
        i0,i1,i2 = edges.searchsorted([wa[0], 1210*zp1, wa[-1]])
        if debug:
            print(i0,i1,i2)
        return edges[i0:i2]
    
    def update_knots(knots, indices, fl, masked):
        """ Calculate the y position of each knot.
    
        Updates `knots` inplace.
    
        Parameters
        ----------
        knots: list of [xpos, ypos, bool] with length N
          bool says whether the knot should kept unchanged.
        indices: list of (i0,i1) index pairs
           The start and end indices into fl and masked of each
           spectrum chunk (xpos of each knot are the chunk centres).
        fl, masked: arrays shape (M,)
           The flux, and boolean arrays showing which pixels are
           masked.
        """
    
        iy, iflag = 1, 2
        for iknot,(i1,i2) in enumerate(indices):
            if knots[iknot][iflag]:
                continue
    
            f0 = fl[i1:i2]
            m0 = masked[i1:i2]
            f1 = f0[~m0]
            knots[iknot][iy] = np.median(f1)
    
    def linear_co(wa, knots):
        """linear interpolation through the spline knots.
    
        Add extra points on either end to give
        a nice slope at the end points."""
        wavc, mfl = list(zip(*knots))[:2]
        extwavc = ([wavc[0] - (wavc[1] - wavc[0])] + list(wavc) +
                   [wavc[-1] + (wavc[-1] - wavc[-2])])
        extmfl = ([mfl[0] - (mfl[1] - mfl[0])] + list(mfl) +
                  [mfl[-1] + (mfl[-1] - mfl[-2])])
        co = np.interp(wa, extwavc, extmfl)
        return co
    
    def Akima_co(wa, knots):
        """Akima interpolation through the spline knots."""
        x,y,_ = zip(*knots)
        spl = interpolate.Akima1DInterpolator(x, y)
        return spl(wa)
    
    def remove_bad_knots(knots, indices, masked, fl, er, debug=False):
        """ Remove knots in chunks without any good pixels. Modifies
        inplace."""
        idelknot = []
        for iknot,(i,j) in enumerate(indices):
            if np.all(masked[i:j]) or np.median(fl[i:j]) <= 2*np.median(er[i:j]):
                if debug:
                    print('Deleting knot', iknot, 'near {:.1f} Angstroms'.format(
                        knots[iknot][0]))
                idelknot.append(iknot)
    
        for i in reversed(idelknot):
            del knots[i]
            del indices[i]
    
    def chisq_chunk(model, fl, er, masked, indices, knots, chithresh=1.5):
        """ Calc chisq per chunk, update knots flags inplace if chisq is
        acceptable. """
        chisq = []
        FLAG = 2
        for iknot,(i1,i2) in enumerate(indices):
            if knots[iknot][FLAG]:
                continue
    
            f0 = fl[i1:i2]
            e0 = er[i1:i2]
            m0 = masked[i1:i2]
            f1 = f0[~m0]
            e1 = e0[~m0]
            mod0 = model[i1:i2]
            mod1 = mod0[~m0]
            resid = (mod1 - f1) / e1
            chisq = np.sum(resid*resid)
            rchisq = chisq / len(f1)
            if rchisq < chithresh:
                #print (good reduced chisq in knot', iknot)
                knots[iknot][FLAG] = True
    
    def prepare_knots(wa, fl, er, edges, ax=None, debug=False):
        """ Make initial knots for the continuum estimation.
    
        Parameters
        ----------
        wa, fl, er : arrays
           Wavelength, flux, error.
        edges : The edges of the wavelength chunks. Splines knots are to be
           places at the centre of these chunks.
        ax : Matplotlib Axes
           If not None, use to plot debugging info.
    
        Returns
        -------
        knots, indices, masked
          * knots: A list of [x, y, flag] lists giving the x and y position
            of each knot.
          * indices: A list of tuples (i,j) giving the start and end index
            of each chunk.
          * masked: An array the same shape as wa.
        """
        indices = wa.searchsorted(edges)
        indices = [(i0,i1) for i0,i1 in zip(indices[:-1],indices[1:])]
        wavc = [0.5*(w1 + w2) for w1,w2 in zip(edges[:-1],edges[1:])]
    
        knots = [[wavc[i], 0, False] for i in range(len(wavc))]
    
        masked = np.zeros(len(wa), bool)
        masked[~(er > 0)] = True
    
        # remove bad knots
        remove_bad_knots(knots, indices, masked, fl, er, debug=debug)
    
        if ax is not None:
            yedge = np.interp(edges, wa, fl)
            ax.vlines(edges, 0, yedge + 100, color='c', zorder=10)
    
        # set the knot flux values
        update_knots(knots, indices, fl, masked)
    
        if ax is not None:
            x,y = list(zip(*knots))[:2]
            ax.plot(x, y, 'o', mfc='none', mec='c', ms=10, mew=1, zorder=10)
    
        return knots, indices, masked
    
    
    def unmask(masked, indices, wa, fl, er, minpix=3):
        """ Forces each chunk to use at least minpix pixels.
    
        Sometimes all pixels can become masked in a chunk. We don't want
        this! This forces there to be at least minpix pixels used in each
        chunk.
         """
        for iknot,(i,j) in enumerate(indices):
            #print(iknot, wa[i], wa[j], (~masked[i:j]).sum())
            if np.sum(~masked[i:j]) < minpix:
                #print('unmasking pixels')
                # need to unmask minpix
                f0 = fl[i:j]
                e0 = er[i:j]
                ind = np.arange(i,j)
                f1 = f0[e0 > 0]
                isort = np.argsort(f1)
                ind1 = ind[e0 > 0][isort[-minpix:]]
                #    print(wa[i], wa[j])
                #    print(wa[ind1])
                masked[ind1] = False
    

    def estimate_continuum(s, knots, indices, masked, ax=None, maxiter=1000,
                           nsig=1.5, debug=False):
        """ Iterate to estimate the continuum.
        """
        count = 0
        while True:
            if debug:
                print('iteration', count)
            update_knots(knots, indices, s.fl, masked)
            model = linear_co(s.wa, knots)
            model_a = Akima_co(s.wa, knots)
            chisq_chunk(model_a, s.fl, s.er, masked,
                        indices, knots, chithresh=1)
            flags = list(zip(*knots))[-1]
            if np.all(flags):
                if debug:
                    print('All regions have satisfactory fit, stopping')
                break
            # remove outliers
            c0 = ~masked
            resid = (model - s.fl) / s.er
            oldmasked = masked.copy()
            masked[(resid > nsig) & ~masked] = True
            unmask(masked, indices, s.wa, s.fl, s.er)
            if np.all(oldmasked == masked):
                if debug:
                    print('No further points masked, stopping')
                break
            if count > maxiter:
                raise RuntimeError('Exceeded maximum iterations')
    
            count +=1

        co = Akima_co(s.wa, knots)
        c0 = co <= 0
        co[c0] = 0
    
        if ax is not None:
            ax.plot(s.wa, linear_co(s.wa, knots), color='0.7', lw=2)
            ax.plot(s.wa, co, 'k', lw=2, zorder=10)
            x,y = list(zip(*knots))[:2]
            ax.plot(x, y, 'o', mfc='none', mec='k', ms=10, mew=1, zorder=10)
    
        return co

    ### Here starts the actual fitting
    ## Pull uncertainty from spectrum
    ## TODO this is very hacky right now
    if not hasattr(spec, "uncertainty"):
        logging.info("No uncertainty, assuming all are equal (continuum will probably fail)")
        error = np.ones(len(spec.wavelength.value))
    else:
        if isinstance(spec.uncertainty, StdDevUncertainty):
            error = spec.uncertainty.array
        else:
            raise ValueError("Could not understand uncertainty type: {}".format(
                    spec.uncertainty))
            
    s = np.rec.fromarrays([spec.wavelength.value,
                           spec.flux.value,
                           error], names=["wa","fl","er"])

    if edges is not None:
        edges = list(edges)
    elif kind.upper() == 'QSO':
        if 'redshift' in kwargs:
            z = kwargs['redshift']
        elif 'redshift' in spec.meta:
            z = spec.meta['redshift']
        else:
            raise RuntimeError(
                "I need the emission redshift for kind='qso'; please\
                provide redshift using `redshift` keyword.")

        divmult = kwargs.get('divmult', 2)
        forest_divmult = kwargs.get('forest_divmult', 2)
        edges = make_chunks_qso(
            s.wa, z, debug=debug, divmult=divmult,
            forest_divmult=forest_divmult)
    
    
    if ax is not None:
        ax.plot(s.wa, s.fl, '-', color='0.4', drawstyle='steps-mid')
        ax.plot(s.wa, s.er, 'g')

    knots, indices, masked = prepare_knots(s.wa, s.fl, s.er, edges,
                                           ax=ax, debug=debug)

    # Note this modifies knots and masked inplace
    co = estimate_continuum(s, knots, indices, masked, ax=ax, debug=debug)

    if ax is not None:
        ax.plot(s.wa[~masked], s.fl[~masked], '.y')
        ymax = np.percentile(s.fl[~np.isnan(s.fl)],  95)
        ax.set_ylim(-0.02*ymax, 1.1*ymax)

    return co, [k[:2] for k in knots]
