import logging
from collections import MutableSequence

import astropy.units as u
import numpy as np
from astropy.nddata import NDIOMixin, NDUncertainty, NDDataRef

from .spectrum1d import Spectrum1D

__all__ = ['SpectrumCollection']


class SpectrumCollection(NDIOMixin):
    """
    A container class for :class:`~specutils.Spectrum1D` objects. This allows
    for operations to be performed over a set of spectrum objects. This class
    behaves similarly to a :class:`~specutils.Spectrum1D`, providing properties
    off the input spectra as n-dimensional :class:`~astropy.units.Quantity`.

    Notes
    -----
    Items in this collection must currently be the same shape. Items are not
    automatically resampled onto a shared grid. Currently, users should resample
    their spectra prior to creating a :class:`~specutils.SpectrumCollection`
    object.

    Parameters
    ----------
    spectra : list, ndarray
        A list of :class:`~specutils.Spectrum1D` objects to be held in the
        collection.
    """
    def __init__(self, flux, spectral_axis, wcs=None, uncertainty=None, 
                 mask=None, meta=None):
        # Ensure that the input values are the same shape
        if not (flux.shape == spectral_axis.shape == wcs.shape == 
                mask.shape == meta.shape):
            raise Exception("Shape of all elements must be the same.")

        # Check for quantity
        if not isinstance(flux, u.Quantity):
            raise u.UnitsError("Flux must be a `Quantity`.")

        if not isinstance(spectral_axis, u.Quantity):
            raise u.UnitsError("Spectral axis must be a `Quantity`.")

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

    @staticmethod
    def from_spectra(spectra):
        # Enforce that the shape of each item must be the same
        if not all((x.shape == spectra[0].shape for x in spectra)):
            raise Exception("Shape of all elements must be the same.")

        # Compose multi-dimensional ndarrays for each property
        flux = np.vstack(
            [spec.flux for spec in spectra]) * spectra[0].flux.unit
        spectral_axis = np.vstack(
            [spec.spectral_axis for spec in spectra]) * spectra[0].spectral_axis.unit
        uncertainty = spectra[0].uncertainty.__class__(
            np.vstack([spec.uncertainty.array for spec in spectra]))
        mask = np.vstack([spec.mask for spec in spectra])

        # Store the wcs and meta as lists
        wcs = [spec.wcs for spec in spectra]
        meta = [spec.meta for spec in spectra]

        return SpectrumCollection(flux, spectral_axis, uncertainty=uncertainty,
                                  wcs=wcs, mask=mask, meta=meta)

    @property
    def flux(self):
        return self._flux

    @property
    def spectral_axis(self):
        return self._spectral_axis

    @property
    def frequency(self):
        return self.spectral_axis.to(u.GHz, u.spectral())

    @property
    def wavelength(self):
        return self.spectral_axis.to(u.AA, u.spectral())

    @property
    def energy(self):
        return self.spectral_axis.to(u.eV, u.spectral())

    @property
    def shape(self):
        return self.flux.shape

    @property
    def ndim(self):
        return self.flux.ndim

    @property
    def uncertainty(self):
        return self._uncertainty


class MutableSpectrumCollection(MutableSequence, NDIOMixin):
    """
    A container class for :class:`~specutils.Spectrum1D` objects. This allows
    for operations to be performed over a set of spectrum objects. This class
    behaves as if it were a :class:`~specutils.Spectrum1D` and supports
    pythonic collection operations like slicing, deleting, and inserting.

    Note
    ----
    Items in this collection must currently be the same shape. Items are not
    automatically resampled onto a shared grid. Currently, users should resample
    their spectra prior to creating a :class:`~specutils.SpectrumCollection`
    object.

    Parameters
    ----------
    items : list, ndarray
        A list of :class:`~specutils.Spectrum1D` objects to be held in the
        collection
    """
    def __init__(self, items):
        # Enforce that the shape of each item must be the same
        if not all((x.shape == items[0].shape for x in items)):
            raise Exception("Shape of all elements must be the same.")

        self._items = items

    def __setitem__(self, index, value):
        self._items[index] = value

    def __delitem__(self, index):
        del self._items[index]

    def insert(self, index, value):
        self._items.insert(index, value)

    def __getitem__(self, index):
        return self._items[index]

    def __len__(self):
        return len(self._items)

    def _get_property(self, name):
        val = [getattr(x, name) for x in self._items]

        if name == 'uncertainty':
            val = [x.array if val is not None else 0 for x in val]

            logging.info(
                "`SpectrumCollection` assumes that all "
                "spectra have the same uncertainty type.")

            if self._items[0].uncertainty is not None:
                val = self._items[0].uncertainty.__class__(val)

        elif hasattr(val[0], 'unit'):
            val = np.vstack(np.array(val)) * val[0].unit

        return val

    def __getattr__(self, name):
        """
        This is a proxy function that forces `SpectrumCollection` to behave
        like a `Spectrum1D` object without having to manually provide
        transient attributes that match the `Spectrum1D` API.
        """
        # TODO: currently, this method assumes that all uncertainties share
        # the same uncertainty type.
        if hasattr(Spectrum1D, name):
            # First, check to see if this property is cached
            val = self._get_property(name)

            return val

        return object.__getattribute__(self, name)

    def __repr__(self):
        return """<SpectrumCollection(size={})>""".format(len(self))
