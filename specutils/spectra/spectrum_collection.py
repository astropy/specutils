import warnings

import astropy.units as u
import numpy as np
from astropy.nddata import NDUncertainty, StdDevUncertainty
from astropy.coordinates import SpectralCoord

from .spectrum1d import Spectrum1D
from astropy.nddata import NDIOMixin

__all__ = ['SpectrumCollection']


class SpectrumCollection(NDIOMixin):
    """
    A class to represent a heterogeneous set of spectra that are the same length
    but have different spectral axes. Spectra that meet this requirement can be
    stored as multidimensional arrays, and thus can have operations performed
    on them faster than if they are treated as individual
    :class:`~specutils.Spectrum1D` objects.

    For multidimensional spectra that all have the *same* spectral axis, use a
    :class:`~specutils.Spectrum1D` with dimension greater than 1.  For a
    collection of spectra that have different shapes, use
    :class:`~specutils.SpectrumList`. For more on this topic, see
    :ref:`specutils-representation-overview`.

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
    wcs : wcs object or None
        A wcs object (if available) for the collection of spectra.  The object
        must follow standard indexing rules to get a sub-wcs if it is provided.
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
            spectral_axis = SpectralCoord(spectral_axis)

            # Ensure that the input values are the same shape
            if not (flux.shape == spectral_axis.shape):
                raise ValueError("Shape of all data elements must be the same.")

        if uncertainty is not None and uncertainty.array.shape != flux.shape:
            raise ValueError("Uncertainty must be the same shape as flux and "
                            "spectral axis.")

        if mask is not None and mask.shape != flux.shape:
            raise ValueError("Mask must be the same shape as flux and "
                            "spectral axis.")

        # Convert uncertainties to standard deviations if not already defined
        # to be of some type
        if uncertainty is not None and not isinstance(uncertainty, NDUncertainty):
            # If the uncertainties are not provided a unit, raise a warning
            # and use the flux units
            if not isinstance(uncertainty, u.Quantity):
                warnings.warn("No unit associated with uncertainty, assuming"
                              f"flux units of '{flux.unit}'.")
                uncertainty = u.Quantity(uncertainty, unit=flux.unit)

            uncertainty = StdDevUncertainty(uncertainty)

        self._flux = flux
        self._spectral_axis = spectral_axis
        self._wcs = wcs
        self._uncertainty = uncertainty
        self._mask = mask
        self._meta = meta

    def __getitem__(self, key):
        flux = self.flux[key]
        if flux.ndim != 1:
            raise ValueError("Currently only 1D data structures may be "
                             "returned from slice operations.")
        spectral_axis = self.spectral_axis[key]
        uncertainty = None if self.uncertainty is None else self.uncertainty[key]
        wcs = None if self.wcs is None else self.wcs[key]
        mask = None if self.mask is None else self.mask[key]
        if self.meta is None:
            meta = None
        else:
            try:
                meta = self.meta[key]
            except KeyError:
                meta = self.meta

        return Spectrum1D(flux=flux, spectral_axis=spectral_axis,
                          uncertainty=uncertainty, wcs=wcs, mask=mask,
                          meta=meta)

    @classmethod
    def from_spectra(cls, spectra):
        """
        Create a spectrum collection from a set of individual
        :class:`specutils.Spectrum1D` objects.

        Parameters
        ----------
        spectra : list, ndarray
            A list of :class:`~specutils.Spectrum1D` objects to be held in the
            collection. Currently the spectral_axis parameters (e.g. observer,
            radial_velocity) must be the same for each spectrum.
        """
        # Enforce that the shape of each item must be the same
        if not all((x.shape == spectra[0].shape for x in spectra)):
            raise ValueError("Shape of all elements must be the same.")

        # Compose multi-dimensional ndarrays for each property
        flux = u.Quantity([spec.flux for spec in spectra])

        # Check that the spectral parameters are the same for each input
        # spectral_axis and create the multi-dimensional SpectralCoord
        sa = [x.spectral_axis for x in spectra]
        if (not all(x.radial_velocity == sa[0].radial_velocity for x in sa) or
                not all(x.target == sa[0].target for x in sa) or
                not all(x.observer == sa[0].observer for x in sa) or
                not all(x.doppler_convention == sa[0].doppler_convention for
                        x in sa) or
                not all(x.doppler_rest == sa[0].doppler_rest for x in sa)):
            raise ValueError("All input spectral_axis SpectralCoord "
                             "objects must have the same parameters.")
        spectral_axis = SpectralCoord(sa,
                            radial_velocity=sa[0].radial_velocity,
                            doppler_rest=sa[0].doppler_rest,
                            doppler_convention=sa[0].doppler_convention,
                            observer=sa[0].observer,
                            target=sa[0].target)

        # Check that either all spectra have associated uncertainties, or that
        # none of them do. If only some do, log an error and ignore the
        # uncertainties.
        if (not all((x.uncertainty is None for x in spectra)) and
                any((x.uncertainty is not None for x in spectra)) and
                all((x.uncertainty.uncertainty_type ==
                     spectra[0].uncertainty.uncertainty_type
                     for x in spectra))):
            quncs = u.Quantity([spec.uncertainty.quantity for spec in spectra])
            uncertainty = spectra[0].uncertainty.__class__(quncs)
        else:
            uncertainty = None
            warnings.warn("Not all spectra have associated uncertainties of "
                          "the same type, skipping uncertainties.")

        # Check that either all spectra have associated masks, or that
        # none of them do. If only some do, log an error and ignore the masks.
        if (not all((x.mask is None for x in spectra)) and
                any((x.mask is not None for x in spectra))):
            mask = np.array([spec.mask for spec in spectra])
        else:
            mask = None
            warnings.warn("Not all spectra have associated masks, "
                          "skipping masks.")

        # Store the wcs and meta as lists
        wcs = [spec.wcs for spec in spectra]
        meta = [spec.meta for spec in spectra]

        return cls(flux=flux, spectral_axis=spectral_axis,
                   uncertainty=uncertainty, wcs=wcs, mask=mask, meta=meta)

    @property
    def flux(self):
        """The flux in the spectrum as a `~astropy.units.Quantity`."""
        return self._flux

    @property
    def spectral_axis(self):
        """The spectral axes as a `~astropy.units.Quantity`."""
        return self._spectral_axis

    @property
    def frequency(self):
        """
        The spectral axis as a frequency `~astropy.units.Quantity` (in GHz).
        """
        return self.spectral_axis.to(u.GHz, u.spectral())

    @property
    def wavelength(self):
        """
        The spectral axis as a wavelength `~astropy.units.Quantity` (in
        Angstroms).
        """
        return self.spectral_axis.to(u.AA, u.spectral())

    @property
    def energy(self):
        """
        The spectral axis as an energy `~astropy.units.Quantity` (in eV).
        """
        return self.spectral_axis.to(u.eV, u.spectral())

    @property
    def wcs(self):
        """The WCS's as an object array"""
        return self._wcs

    @property
    def uncertainty(self):
        """The uncertainty in the spectrum as a `~astropy.units.Quantity`."""
        return self._uncertainty

    @property
    def mask(self):
        """The mask array for the spectrum."""
        return self._mask

    @property
    def meta(self):
        """A dictionary of metadata for theis spectrum collection, or `None`."""
        return self._meta

    @property
    def shape(self):
        """
        The shape of the collection. This is *not* the same as
        the shape of the flux et al., because the trailing (spectral)
        dimension is not included here.
        """
        return self.flux.shape[:-1]

    def __len__(self):
        return self.shape[0]

    @property
    def ndim(self):
        """
        The dimensionality of the collection. This is *not* the same as
        the dimensionality of the flux et al., because the trailing (spectral)
        dimension is not included here.
        """
        return self.flux.ndim - 1

    @property
    def nspectral(self):
        """
        The length of the spectral dimension.
        """
        return self.flux.shape[-1]

    def __repr__(self):
        return """SpectrumCollection(ndim={}, shape={})
    Flux units:          {}
    Spectral axis units: {}
    Uncertainty type:    {}""".format(
            self.ndim, self.shape, self.flux.unit, self.spectral_axis.unit,
            self.uncertainty.uncertainty_type if self.uncertainty is not None else None)
