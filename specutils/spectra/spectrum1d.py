from __future__ import division

import numpy as np
from astropy import units as u
from astropy.nddata import NDDataRef
from astropy.utils.decorators import lazyproperty

from ..wcs import WCSWrapper
from .spectrum_mixin import OneDSpectrumMixin

__all__ = ['Spectrum1D']


class Spectrum1D(OneDSpectrumMixin, NDDataRef):
    """
    Spectrum container for 1D spectral data.
    """

    def __init__(self, flux, spectral_axis=None, wcs=None, unit=None,
                 spectral_axis_unit=None, *args, **kwargs):
        # Attempt to parse the WCS. If not WCS object is given, try instead to
        # parse a given wavelength array. This is put into a GWCS object to
        # then be used behind-the-scenes for all specutils operations.
        if wcs is not None:
            wcs = WCSWrapper(wcs)
        elif spectral_axis is not None:
            spectral_axis = u.Quantity(spectral_axis,
                                       unit=spectral_axis_unit)

            wcs = WCSWrapper.from_array(spectral_axis)
        else:
            # If not wcs and not spectral axis has been given, raise an error
            raise LookupError("No WCS object or spectral axis information has "
                              "been given. Please provide one.")

        if not isinstance(flux, u.Quantity):
            flux = u.Quantity(flux, unit=unit or "Jy")

        self._velocity_convention = None

        # Currently, only a fits wcs object stores the rest wavelength or
        # frequency information in the wcs object. In any other case, the user
        # will be given an warning to provide these explicitly in this object.
        if wcs.rest_frequency != 0:
            self._rest_value = wcs.rest_frequency * u.Hz
        elif wcs.rest_wavelength != 0:
            self._rest_value = wcs.rest_wavelength * u.AA

        super(Spectrum1D, self).__init__(data=flux.value, unit=flux.unit,
                                         wcs=wcs, *args, **kwargs)

    @property
    def frequency(self):
        return self.spectral_axis.to(u.GHz, u.spectral())

    @property
    def wavelength(self):
        return self.spectral_axis.to(u.AA, u.spectral())

    @property
    def energy(self):
        return self.spectral_axis.to(u.eV, u.spectral())

    @lazyproperty
    def bin_edges(self):
        return self.wcs.bin_edges()

    @velocity_convention.setter
    def velocity_convention(self, value):
        if value not in ('relativistic', 'optical', 'radio'):
            raise ValueError("The allowed velocity conventions are 'optical' "
                             "(linear with respect to wavelength), 'radio' "
                             "(linear with respect to frequency), and "
                             "'relativistic'.")
        self._velocity_convention = value

    @rest_value.setter
    def rest_value(self, value):
        if not hasattr(value, 'unit') or not value.unit.is_equivalent(u.Hz, u.spectral()):
            raise ValueError(
                "Rest value must be energy/wavelength/frequency equivalent.")

        self._rest_value = value

    def spectral_resolution(self, true_dispersion, delta_dispersion, axis=-1):
        """Evaluate the probability distribution of the spectral resolution.

        For example, to tabulate a binned resolution function at 6000A
        covering +/-10A in 0.2A steps::

        R = spectrum1d.spectral_resolution(
        ... 6000 * u.Angstrom, np.linspace(-10, 10, 51) * u.Angstrom)
        assert R.shape == (50,)
        assert np.allclose(R.sum(), 1.)

        To build a sparse resolution matrix for true wavelengths 4000-8000A
        in 0.1A steps::

        R = spectrum1d.spectral_resolution(
        ... np.linspace(4000, 8000, 40001)[:, np.newaxis] * u.Angstrom,
        ... np.linspace(-10, +10, 201) * u.Angstrom)
        assert R.shape == (40000, 200)
        assert np.allclose(R.sum(axis=1), 1.)

        Parameters
        ----------
        true_dispersion : ~`astropy.units.Quantity`
            True value(s) of dispersion for which the resolution should be
            evaluated.
        delta_dispersion : ~`astropy.units.Quantity`
            Array of (observed - true) dispersion bin edges to integrate the
            resolution probability density over.
        axis : int
            Which axis of ``delta_dispersion`` contains the strictly increasing
            dispersion values to interpret as bin edges.  The dimension of
            ``delta_dispersion`` along this axis must be at least two.

        Returns
        -------
        numpy array
            Array of dimensionless probabilities calculated as the integral of
            P(observed | true) over each bin in (observed - true). The output
            shape is the result of broadcasting the input shapes.
        """
        pass
