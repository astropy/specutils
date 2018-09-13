import logging

import astropy.units as u
import numpy as np
from astropy.nddata import NDUncertainty, StdDevUncertainty

from .spectrum1d import Spectrum1D

__all__ = ['SpectrumCollection']


class SpectrumCollection:
    """
    A class to represent a heterogeneous set of spectra that are the same length
    but have different spectral axes. Spectra that meet this requirement can be
    stored as multidimensional arrays, and thus can have operations performed
    on them faster than if they are treated as individual
    :class:`~specutils.Spectrum1D` objects.

    The attributes on this class uses the same names and conventions as
    :class:`~specutils.Spectrum1D`, allowing some operations to work the same.
    Where this does not work, the user can use standard indexing notation to
    access individual :class:`~specutils.Spectrum1D` objects.

    Parameters
    ----------
    flux : :class:`astropy.units.Quantity`
        The flux data. The trailing dimension should be the spectral dimension.
    spectral_axis : :class:`astropy.units.Quantity`
        The spectral axes of the spectra (e.g., wavelength).  Must match the
        dimensionality of ``flux``.
    wcs : list or None
        A list of the input WCS associated with each set of spectrum of the
        collection (if needed).
    uncertainty : :class:`astropy.nddata.NDUncertainty` or ndarray
        The uncertainties associated with each spectrum of the collection. In
        the case that only an n-dimensional quantity or ndaray is provided,
        the uncertainties are assumed to be standard deviations. Must match the
        dimensionality of ``flux``.
    mask : ndarray or None
        The n-dimensional mask information associated with each spectrum. If
        present, must match the dimensionality of ``flux``.
    meta : list
        The list of dictionaries containing meta data to be associated with
        each spectrum in the collection.
    """
    def __init__(self, flux, spectral_axis=None, wcs=None, uncertainty=None,
                 mask=None, meta=None):
        # Check for quantity
        if not isinstance(flux, u.Quantity):
            raise u.UnitsError("Flux must be a `Quantity`.")

        if spectral_axis is not None:
            if not isinstance(spectral_axis, u.Quantity):
                raise u.UnitsError("Spectral axis must be a `Quantity`.")

            # Ensure that the input values are the same shape
            if not (flux.shape == spectral_axis.shape):
                raise ValueError("Shape of all data elements must be the same.")

        if uncertainty is not None and uncertainty.array.shape != flux.shape:
            raise ValueError("Uncertainty must be the same shape as flux and "
                            "spectral axis.")

        if mask is not None and mask.shape != flux.shape:
            raise ValueError("Mask must be the same shape as flux and "
                            "spectral axis.")

        if wcs is not None:
            wcs = np.array(wcs, copy=False, dtype=object)
            if wcs.shape != flux.shape:
                raise ValueError("All spectra must be associated with a single WCS.")

        # Convert uncertainties to standard deviations if not already defined
        # to be of some type
        if not isinstance(uncertainty, NDUncertainty):
            # If the uncertainties are not provided a unit, raise a warning
            # and use the flux units
            if not isinstance(uncertainty, u.Quantity):
                logging.warning("No unit associated with uncertainty, assuming"
                                "flux units of '{}'.".format(flux.unit))
                uncertainty = u.Quantity(uncertainty, unit=flux.unit)

            uncertainty = StdDevUncertainty(uncertainty)

        self._flux = flux
        self._spectral_axis = spectral_axis
        self._wcs = wcs
        self._uncertainty = uncertainty
        self._mask = mask
        self._meta = meta

    def __getitem__(self, key):
        if isinstance(key, int) and self.ndim > 2:
            raise ValueError("Currently only 1D data structures may be "
                             "returned from slice operations.")

        return Spectrum1D(flux=self.flux[key],
                          spectral_axis=self.spectral_axis[key],
                          uncertainty=self.uncertainty[key],
                          wcs=self.wcs[key],
                          mask=self.mask[key],
                          meta=self.meta[key])

    @classmethod
    def from_spectra(cls, spectra):
        """
        Create a spectrum collection from a set of individual
        :class:`specutils.Spectrum1D` objects.

        Parameters
        ----------
        spectra : list, ndarray
            A list of :class:`~specutils.Spectrum1D` objects to be held in the
            collection.
        """
        # Enforce that the shape of each item must be the same
        if not all((x.shape == spectra[0].shape for x in spectra)):
            raise ValueError("Shape of all elements must be the same.")

        # Compose multi-dimensional ndarrays for each property
        flux = np.vstack(
            [spec.flux.value for spec in spectra]) * spectra[0].flux.unit
        spectral_axis = np.vstack(
            [spec.spectral_axis.value for spec in spectra]) * spectra[0].spectral_axis.unit

        # Check that either all spectra have associated uncertainties, or that
        # none of them do. If only some do, log an error and ignore the
        # uncertainties.
        if not all((x.uncertainty is None for x in spectra)) and \
            any((x.uncertainty is not None for x in spectra)):
            uncertainty = spectra[0].uncertainty.__class__(
                np.vstack([spec.uncertainty.array for spec in spectra]))
        else:
            uncertainty = None

            logging.warning("Not all spectra have associated uncertainties, "
                            "skipping uncertainties.")

        # Check that either all spectra have associated masks, or that
        # none of them do. If only some do, log an error and ignore the masks.
        if not all((x.mask is None for x in spectra)) and \
            any((x.mask is not None for x in spectra)):
            mask = np.vstack([spec.mask for spec in spectra])
        else:
            mask = None

            logging.warning("Not all spectra have associated masks, "
                            "skipping masks.")

        # Store the wcs and meta as lists
        wcs = [spec.wcs for spec in spectra]
        meta = [spec.meta for spec in spectra]

        return cls(flux=flux, spectral_axis=spectral_axis,
                   uncertainty=uncertainty, wcs=wcs, mask=mask, meta=meta)

    @property
    def flux(self):
        """Return the n-dimensional flux quantity object."""
        return self._flux

    @property
    def spectral_axis(self):
        """Return the n-dimensional spectral axis quantity object."""
        return self._spectral_axis

    @property
    def frequency(self):
        """Converts and returns the spectral axis in frequency space."""
        return self.spectral_axis.to(u.GHz, u.spectral())

    @property
    def wavelength(self):
        """Converts and returns the spectral axis in wavelength space."""
        return self.spectral_axis.to(u.AA, u.spectral())

    @property
    def energy(self):
        """Converts and returns the spectral axis in energy space."""
        return self.spectral_axis.to(u.eV, u.spectral())

    @property
    def wcs(self):
        """Return the list of wcs objects."""
        return self._wcs

    @property
    def uncertainty(self):
        """Return the n-dimensional uncertainty object."""
        return self._uncertainty

    @property
    def mask(self):
        """Return the n-dimesional mask object, or `None`."""
        return self._mask

    @property
    def meta(self):
        """Return the list of meta dictionary objects, or `None`."""
        return self._meta

    @property
    def shape(self):
        """Get the shape of the collection."""
        return self.flux.shape

    @property
    def ndim(self):
        """Get the dimensionality of the collection."""
        return self.flux.ndim

    def __repr__(self):
        return """SpectrumCollection(ndim={}, shape={})
    Flux units:          {}
    Spectral axis units: {}
    Uncertainty type:    {}""".format(
            self.ndim, self.shape, self.flux.unit, self.spectral_axis.unit,
            self.uncertainty.uncertainty_type if self.uncertainty is not None else None)
