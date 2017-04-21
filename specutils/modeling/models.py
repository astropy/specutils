from __future__ import absolute_import, division
from astropy.extern import six

from astropy.modeling import Fittable1DModel, Parameter, models
from astropy.modeling.core import _ModelMeta
from astropy.modeling.models import Linear1D
import numpy as np

from .profiles import TauProfile


class Voigt1D(Fittable1DModel):
    """
    Implements a Voigt profile (convolution of Cauchy-Lorentz and Gaussian
    distribution).
    """
    lambda_0 = Parameter(min=0)
    f_value = Parameter(min=1e-4, max=2.0, fixed=True)
    gamma = Parameter(fixed=True)
    v_doppler = Parameter(default=1e5)
    column_density = Parameter(default=13)
    delta_v = Parameter(default=0, fixed=True)
    delta_lambda = Parameter(default=0, fixed=True)

    def evaluate(self, x, lambda_0, f_value, gamma, v_doppler, column_density,
                 delta_v, delta_lambda):
        profile = TauProfile(x, lambda_0=lambda_0, f_value=f_value,
                             gamma=gamma, v_doppler=v_doppler,
                             column_density=column_density,
                             n_lambda=x.size,
                             delta_v=delta_v, delta_lambda=delta_lambda)

        flux = np.exp(-profile.optical_depth) - 1.0

        return flux

    @staticmethod
    def fit_deriv(x, x_0, b, gamma, f):
        return [0, 0, 0, 0]


class AbsorptionMeta(type):
    """
    Meta class for allowing arbitrary numbers of absorption line features in
    spectral models.
    """
    def __call__(cls, lines=None, continuum=None):
        """
        See the `Absorption1D` initializer for details.
        """
        mod_list = []

        if continuum is not None:
            if issubclass(continuum.__class__, Fittable1DModel):
                mod_list.append(continuum)
            elif isinstance(continuum, six.string_types):
                continuum = getattr(models, continuum, 'Linear1D')
                mod_list.append(continuum)
            else:
                raise AttributeError("Unknown continuum type {}.".format(
                    type(continuum)))
        else:
            continuum = Linear1D(slope=0, intercept=1)
            mod_list.append(continuum)

        if lines is not None:
            if isinstance(lines, list):
                mod_list += lines
            elif issubclass(lines, Fittable1DModel):
                mod_list.append(lines)

        abs_mod = np.sum(mod_list)

        def call(self, dispersion, uncertainty=None, dispersion_unit=None,
                 *args, **kwargs):
            from ..spectra import Spectrum1D

            lines = abs_mod._submodels[1:]  # noqa

            flux = super(abs_mod.__class__, self).__call__(dispersion,
                                                           *args, **kwargs)
            spectrum = Spectrum1D(flux, dispersion=dispersion,
                                  dispersion_unit=dispersion_unit)

            return spectrum

        mod = type('Absorption1D', (abs_mod.__class__, ), {'__call__': call})

        return mod()


@six.add_metaclass(type('CombinedMeta', (_ModelMeta, AbsorptionMeta), {}))
class Absorption1D(Fittable1DModel):
    """
    One dimensional spectral model that allows the addition of an arbitrary
    number of absorption features.
    """
    def __init__(self, lines=None, continuum=None, *args, **kwargs):
        """
        Custom fittable model representing a spectrum object.

        Parameters
        ----------
        lines : list
            List containing spectral line features as instances of the
            :class:`~spectacle.core.lines.Line` class.

        continuum : :class:`~astropy.modeling.Fittable1DModel` or str
            The `continuum` argument can either be a model instance, or a
            string representing a model type (see `Astropy models list
            <http://docs.astropy.org/en/stable/modeling/index.html#
            module-astropy.modeling.functional_models>`_ for options.

            .. note:: Continuum classes instantiating by string reference will
                      be initialized with default parameters.
        """
        super(Absorption1D, self).__init__(*args, **kwargs)

    def get_range_mask(self, dispersion, name=None):
        """
        Returns a mask for the model spectrum that indicates where the values
        are discernible from the continuum. I.e. where the absorption data
        is contained.

        If a `x_0` value is provided, only the related Voigt profile will be
        considered, otherwise, the entire model is considered.

        Parameters
        ----------
        dispersion : array-like
            The x values for which the model spectrum will be calculated.
        name : str, optional
            The line name used to grab a specific Voigt profile.

        Returns
        -------
        array-like
            Boolean array indicating indices of interest.
        """
        profile = self.model if name is None else self.get_profile(name)
        vdisp = profile(dispersion)
        cont = np.zeros(dispersion.shape)

        return ~np.isclose(vdisp, cont, rtol=1e-2, atol=1e-5)

    def get_profile(self, name):
        """
        Retrieve the particular :class:`spectacle.modeling.models.Voigt1D`
        model for a given absorption feature.

        Parameters
        ----------
        name : str
            The unique name of the model.

            .. note:: Name uniqueness is not enforced. If there is more than
                      one line with the same name, only the first will be
                      returned.

        Returns
        -------
        :class:`spectacle.modeling.models.Voigt1D`
            The Voigt model for the particular line.
        """
        return next((sm for sm in self._submodels if sm.name == name), None)

    def get_fwhm(self, line_name):
        """
        Return the full width at half max for a given absorption line feature.

        Parameters
        ----------
        line_name : str
            The name of the absorption line feature.

        Returns
        -------
        float
            The calculated full width at half max.
        """
        profile = self.get_profile(line_name)

        return profile.fwhm
