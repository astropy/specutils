from __future__ import division

import astropy.units as u
from astropy import constants as c
import numpy as np


class TauProfile(object):
    """
    Create an optical depth vs. wavelength profile for an
    absorption line using a voigt profile. This follows the paradigm of
    the :func:`~yt.analysis_modules.absorption_spectrum.absorption_line`
    profile generator.

    Parameters
    ----------
    lambda_0 : float
       Central wavelength in Angstroms.
    f_value : float
       Absorption line oscillator strength.
    gamma : float
       Absorption line gamma value.
    v_doppler : float
       Doppler b-parameter in cm/s.
    column_density : float
       Column density in cm^-2.
    delta_v : float
       Velocity offset from lambda_0 in cm/s. Default: None (no shift).
    delta_lambda : float
        Wavelength offset in Angstrom. Default: None (no shift).
    lambda_bins : array-like
        Wavelength array for line deposition in Angstroms. If None, one will be
        created using n_lambda and dlambda. Default: None.
    n_lambda : int
        Size of lambda bins to create if lambda_bins is None. Default: 12000.
    dlambda : float
        Lambda bin width in Angstroms if lambda_bins is None. Default: 0.01.
    """
    def __init__(self, x, lambda_0, f_value, gamma, v_doppler, column_density,
                 delta_v=None, delta_lambda=None, lambda_bins=None,
                 n_lambda=12000, dlambda=0.01):
        charge_proton = u.Quantity(4.8032056e-10, 'esu')
        tau_factor = ((np.sqrt(np.pi) * charge_proton ** 2 /
                       (u.M_e.cgs * c.c.cgs))).cgs

        # Make the input parameters quantities so we can keep track
        # of units
        x = x * u.Unit('Angstrom')
        lambda_0 = lambda_0 * u.Unit('Angstrom')
        v_doppler = v_doppler * u.Unit('cm/s')
        column_density = 10 ** column_density * u.Unit('1/cm2')
        delta_v = delta_v * u.Unit('cm/s')
        delta_lambda = delta_lambda * u.Unit('Angstrom')

        lambda_bins = x

        # shift lambda_0 by delta_v
        if delta_v is not None:
            lam1 = lambda_0 * (1 + delta_v / c.c.cgs)
        elif delta_lambda is not None:
            lam1 = lambda_0 + delta_lambda
        else:
            lam1 = lambda_0

        # conversions
        nudop = (v_doppler / lam1).to('Hz')  # doppler width in Hz

        # create wavelength
        if lambda_bins is None:
            lambda_bins = lam1 + \
                          np.arange(n_lambda, dtype=np.float) * dlambda - \
                          n_lambda * dlambda / 2  # wave vector (angstroms)

        # tau_0
        tau_X = tau_factor * column_density * f_value / v_doppler
        tau0 = (tau_X * lambda_0).decompose()

        # dimensionless frequency offset in units of doppler freq
        x = c.c.cgs / v_doppler * (lam1 / lambda_bins - 1.0)
        a = gamma / (4.0 * np.pi * nudop)  # damping parameter
        phi = self.voigt(a, x)  # line profile
        tau_phi = tau0 * phi  # profile scaled with tau0

        self._lambda_bins = lambda_bins.value
        self._tau_phi = tau_phi.decompose().value

    @property
    def optical_depth(self):
        return self._tau_phi

    @classmethod
    def voigt(cls, a, u):
        from scipy import special

        x = np.asarray(u).astype(np.float64)
        y = np.asarray(a).astype(np.float64)

        return special.wofz(x + 1j * y).real
