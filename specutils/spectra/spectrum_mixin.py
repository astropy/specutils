import logging

import numpy as np
from astropy.wcs import WCSSUB_SPECTRAL
from astropy.units import Unit
from astropy.nddata import NDData
from astropy import units as u

# Use this once in specutils
# from ..utils.wcs_utils import determine_ctype_from_vconv, convert_spectral_axis

from spectral_cube.spectral_axis import determine_ctype_from_vconv, convert_spectral_axis

__all__ = ['Spectrum1D']


class SpectrumMixin(object):

    @property
    def _spectral_axis_numpy_index(self):
        return self.data.ndim - 1 - self.wcs.wcs.spec

    @property
    def _spectral_axis_len(self):
        """
        How many elements are in the spectral dimension?
        """
        return self.data.shape[self._spectral_axis_numpy_index]

    @property
    def _data_with_spectral_axis_last(self):
        """
        Returns a view of the data with the spectral axis last
        """
        if self._spectral_axis_numpy_index == self.data.ndim - 1:
            return self.data
        else:
            return self.data.swapaxes(self._spectral_axis_numpy_index, self.data.ndim - 1)

    @property
    def _data_with_spectral_axis_first(self):
        """
        Returns a view of the data with the spectral axis first
        """
        if self._spectral_axis_numpy_index == 0:
            return self.data
        else:
            return self.data.swapaxes(self._spectral_axis_numpy_index, 0)

    @property
    def spectral_wcs(self):
        return self.wcs.sub([WCSSUB_SPECTRAL])

    @property
    def spectral_axis(self):
        """
        Returns a Quantity array with the values of the spectral axis.
        """

        spectral_wcs = self.spectral_wcs

        # Lim: What if I have wavelength arrays and I don't want WCS conversion?
        # Tom: this is beyond the scope of the prototype work, your question is more
        #      how to make a WCS object that contains a wavelngth array. The mixin
        #      is for NDData which assumes a WCS has been created (not necessarily
        #      a *FITS* WCS, just some transformation object).

        if spectral_wcs.naxis == 0:
            raise TypeError('WCS has no spectral axis')

        spectral_axis = spectral_wcs.all_pix2world(np.arange(self._spectral_axis_len), 0)[0]

        # Try to get the dispersion unit information
        try:
            spectral_unit = self.wcs.wcs.cunit[self.wcs.wcs.spec]
        except AttributeError:
            logging.warning("No spectral_axis unit information in WCS.")
            spectral_unit = Unit("")

        spectral_axis = spectral_axis * spectral_unit

        return spectral_axis

    def with_spectral_units(self, unit, velocity_convention=None,
                            rest_value=None):
        """
        Example method that returns a new cube with updated spectral axis units
        """

        # TODO: avoid the two function calls, make a wrapper?

        out_ctype = determine_ctype_from_vconv(self.wcs.wcs.ctype[self.wcs.wcs.spec],
                                               unit, velocity_convention=velocity_convention)

        new_wcs = convert_spectral_axis(self._wcs, unit, out_ctype,
                                        rest_value=rest_value)

        return self.__class__(self, wcs=new_wcs)

    # Example methods follow to demonstrate how methods can be written to be
    # agnostic of the non-spectral dimensions.

    def substract_background(self, background):
        """
        Proof of concept, this subtracts a background spectrum-wise
        """

        data = self._data_with_spectral_axis_last

        if callable(background):
            # create substractable array
            pass
        elif (isinstance(background, np.ndarray) and
              background.shape == data[-1].shape):
            substractable_continuum = background
        else:
            raise ValueError("background needs to be callable or have the same shape as the spectum")

        data[-1] -= substractable_continuum

    def normalize(self):
        """
        Proof of concept, this normalizes each spectral dimension based
        on a trapezoidal integration.
        """

        # this gets a view - if we want normalize to return a new NDData object
        # then we should make _data_with_spectral_axis_first return a copy.
        data = self._data_with_spectral_axis_first

        dx = np.diff(self.spectral_axis)
        dy = 0.5 * (data[:-1] + data[1:])

        norm = np.sum(dx * dy.transpose(), axis=-1).transpose()

        data /= norm

    def spectral_interpolation(self, spectral_value, flux_unit=None):
        """
        Proof of concept, this interpolates along the spectral dimension
        """

        data = self._data_with_spectral_axis_last

        from scipy.interpolate import interp1d

        interp = interp1d(self.spectral_axis.value, data)

        x = spectral_value.to(self.spectral_axis.unit, equivalencies=u.spectral())
        y = interp(x)

        if self.unit is not None:
            y *= self.unit

        if flux_unit is None:  # Lim: Is this acceptable?
            return y
        else:
            return y.to(flux_unit, equivalencies=u.spectral_density(x))
