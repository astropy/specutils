from copy import deepcopy
import logging

import numpy as np
from astropy import units as u
from astropy.utils.decorators import lazyproperty
from astropy.nddata import NDUncertainty, NDIOMixin, NDArithmeticMixin

from .spectral_axis import SpectralAxis
from .spectrum_mixin import OneDSpectrumMixin
from .spectral_region import SpectralRegion
from ..utils.wcs_utils import gwcs_from_array
from astropy.coordinates import SpectralCoord
from ndcube import NDCube

__all__ = ['Spectrum1D']

log = logging.getLogger(__name__)


class Spectrum1D(OneDSpectrumMixin, NDCube, NDIOMixin, NDArithmeticMixin):
    """
    Spectrum container for 1D spectral data.

    Note that "1D" in this case refers to the fact that there is only one
    spectral axis.  `Spectrum1D` can contain "vector 1D spectra" by having the
    ``flux`` have a shape with dimension greater than 1.  The requirement
    is that the last dimension of ``flux`` match the length of the
    ``spectral_axis``.

    For multidimensional spectra that are all the same shape but have different
    spectral axes, use a :class:`~specutils.SpectrumCollection`.  For a
    collection of spectra that have different shapes, use
    :class:`~specutils.SpectrumList`. For more on this topic, see
    :ref:`specutils-representation-overview`.

    Parameters
    ----------
    flux : `~astropy.units.Quantity` or `~astropy.nddata.NDData`-like
        The flux data for this spectrum. This can be a simple `~astropy.units.Quantity`,
        or an existing `~Spectrum1D` or `~ndcube.NDCube` object.
    spectral_axis : `~astropy.units.Quantity` or `~specutils.SpectralAxis`
        Dispersion information with the same shape as the last (or only)
        dimension of flux, or one greater than the last dimension of flux
        if specifying bin edges.
    wcs : `~astropy.wcs.WCS` or `~gwcs.wcs.WCS`
        WCS information object that either has a spectral component or is
        only spectral.
    velocity_convention : {"doppler_relativistic", "doppler_optical", "doppler_radio"}
        Convention used for velocity conversions.
    rest_value : `~astropy.units.Quantity`
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency, wave number). Describes the rest value
        of the spectral axis for use with velocity conversions.
    redshift
        See `redshift` for more information.
    radial_velocity
        See `radial_velocity` for more information.
    bin_specification : str
        Either "edges" or "centers" to indicate whether the `spectral_axis`
        values represent edges of the wavelength bin, or centers of the bin.
    uncertainty : `~astropy.nddata.NDUncertainty`
        Contains uncertainty information along with propagation rules for
        spectrum arithmetic. Can take a unit, but if none is given, will use
        the unit defined in the flux.
    mask : `~numpy.ndarray`-like
        Array where values in the flux to be masked are those that
        ``astype(bool)`` converts to True. (For example, integer arrays are not
        masked where they are 0, and masked for any other value.)
    meta : dict
        Arbitrary container for any user-specific information to be carried
        around with the spectrum container object.
    """
    def __init__(self, flux=None, spectral_axis=None, wcs=None,
                 velocity_convention=None, rest_value=None, redshift=None,
                 radial_velocity=None, bin_specification=None, **kwargs):
        # Check for pre-defined entries in the kwargs dictionary.
        unknown_kwargs = set(kwargs).difference(
            {'data', 'unit', 'uncertainty', 'meta', 'mask', 'copy',
             'extra_coords'})

        if len(unknown_kwargs) > 0:
            raise ValueError("Initializer contains unknown arguments(s): {}."
                             "".format(', '.join(map(str, unknown_kwargs))))

        # If the flux (data) argument is already a Spectrum1D (as it would
        # be for internal arithmetic operations), avoid setup entirely.
        if isinstance(flux, Spectrum1D):
            super().__init__(flux)
            return

        # Handle initializing from NDCube objects
        elif isinstance(flux, NDCube):
            if flux.unit is None:
                raise ValueError("Input NDCube missing unit parameter")

            # Change the flux array from bare ndarray to a Quantity
            q_flux = flux.data << u.Unit(flux.unit)

            self.__init__(flux=q_flux, wcs=flux.wcs, mask=flux.mask,
                          uncertainty=flux.uncertainty)
            return

        # If the mask kwarg is not passed to the constructor, but the flux array
        # contains NaNs, add the NaN locations to the mask.
        if "mask" not in kwargs and flux is not None:
            nan_mask = np.isnan(flux)
            if nan_mask.any():
                if hasattr(self, "mask"):
                    kwargs["mask"] = np.logical_or(nan_mask, self.mask)
                else:
                    kwargs["mask"] = nan_mask.copy()
            del nan_mask

        # Ensure that the flux argument is an astropy quantity
        if flux is not None:
            if not isinstance(flux, u.Quantity):
                raise ValueError("Flux must be a `Quantity` object.")
            elif flux.isscalar:
                flux = u.Quantity([flux])

        # Ensure that only one or neither of these parameters is set
        if redshift is not None and radial_velocity is not None:
            raise ValueError("Cannot set both radial_velocity and redshift at "
                             "the same time.")

        # In cases of slicing, new objects will be initialized with `data`
        # instead of ``flux``. Ensure we grab the `data` argument.
        if flux is None and 'data' in kwargs:
            flux = kwargs.pop('data')

        # Ensure that the unit information codified in the quantity object is
        # the One True Unit.
        kwargs.setdefault('unit', flux.unit if isinstance(flux, u.Quantity)
                                            else kwargs.get('unit'))

        # In the case where the arithmetic operation is being performed with
        # a single float, int, or array object, just go ahead and ignore wcs
        # requirements
        if (not isinstance(flux, u.Quantity) or isinstance(flux, float)
            or isinstance(flux, int)) and np.ndim(flux) == 0:

            super(Spectrum1D, self).__init__(data=flux, wcs=wcs, **kwargs)
            return

        if rest_value is None:
            if hasattr(wcs, 'rest_frequency') and wcs.rest_frequency != 0:
                rest_value = wcs.rest_frequency * u.Hz
            elif hasattr(wcs, 'rest_wavelength') and wcs.rest_wavelength != 0:
                rest_value = wcs.rest_wavelength * u.AA
            elif hasattr(wcs, 'wcs') and hasattr(wcs.wcs, 'restfrq') and wcs.wcs.restfrq > 0:
                rest_value = wcs.wcs.restfrq * u.Hz
            elif hasattr(wcs, 'wcs') and hasattr(wcs.wcs, 'restwav') and wcs.wcs.restwav > 0:
                rest_value = wcs.wcs.restwav * u.m
            else:
                rest_value = None
        else:
            if not isinstance(rest_value, u.Quantity):
                log.info("No unit information provided with rest value. "
                             "Assuming units of spectral axis ('%s').",
                             spectral_axis.unit)
                rest_value = u.Quantity(rest_value, spectral_axis.unit)
            elif not rest_value.unit.is_equivalent(u.AA, equivalencies=u.spectral()):
                raise u.UnitsError("Rest value must be "
                                   "energy/wavelength/frequency equivalent.")

        # If flux and spectral axis are both specified, check that their lengths
        # match or are off by one (implying the spectral axis stores bin edges)
        if flux is not None and spectral_axis is not None:
            if spectral_axis.shape[0] == flux.shape[-1]:
                if bin_specification == "edges":
                    raise ValueError("A spectral axis input as bin edges"
                        "must have length one greater than the flux axis")
                bin_specification = "centers"
            elif spectral_axis.shape[0] == flux.shape[-1]+1:
                if bin_specification == "centers":
                    raise ValueError("A spectral axis input as bin centers"
                        "must be the same length as the flux axis")
                bin_specification = "edges"
            else:
                raise ValueError(
                    "Spectral axis length ({}) must be the same size or one "
                    "greater (if specifying bin edges) than that of the last "
                    "flux axis ({})".format(spectral_axis.shape[0],
                                            flux.shape[-1]))

        # If a WCS is provided, check that the spectral axis is last and reorder
        # the arrays if not
        if wcs is not None and hasattr(wcs, "naxis"):
            if wcs.naxis > 1:
                temp_axes = []
                phys_axes = wcs.world_axis_physical_types
                for i in range(len(phys_axes)):
                    if phys_axes[i] is None:
                        continue
                    if phys_axes[i][0:2] == "em" or phys_axes[i][0:5] == "spect":
                        temp_axes.append(i)
                if len(temp_axes) != 1:
                    raise ValueError("Input WCS must have exactly one axis with "
                                     "spectral units, found {}".format(len(temp_axes)))

                # Due to FITS conventions, a WCS with spectral axis first corresponds
                # to a flux array with spectral axis last.
                if temp_axes[0] != 0:
                    log.warning("Input WCS indicates that the spectral axis is not"
                                " last. Reshaping arrays to put spectral axis last.")
                    wcs = wcs.swapaxes(0, temp_axes[0])
                    if flux is not None:
                        flux = np.swapaxes(flux, len(flux.shape)-temp_axes[0]-1, -1)
                    if "mask" in kwargs:
                        if kwargs["mask"] is not None:
                            kwargs["mask"] = np.swapaxes(kwargs["mask"],
                                                len(kwargs["mask"].shape)-temp_axes[0]-1, -1)
                    if "uncertainty" in kwargs:
                        if kwargs["uncertainty"] is not None:
                            if isinstance(kwargs["uncertainty"], NDUncertainty):
                                # Account for Astropy uncertainty types
                                unc_len = len(kwargs["uncertainty"].array.shape)
                                temp_unc = np.swapaxes(kwargs["uncertainty"].array,
                                                       unc_len-temp_axes[0]-1, -1)
                                if kwargs["uncertainty"].unit is not None:
                                    temp_unc = temp_unc * u.Unit(kwargs["uncertainty"].unit)
                                kwargs["uncertainty"] = type(kwargs["uncertainty"])(temp_unc)
                            else:
                                kwargs["uncertainty"] = np.swapaxes(kwargs["uncertainty"],
                                                        len(kwargs["uncertainty"].shape) -
                                                        temp_axes[0]-1, -1)


        # Attempt to parse the spectral axis. If none is given, try instead to
        # parse a given wcs. This is put into a GWCS object to
        # then be used behind-the-scenes for all specutils operations.
        if spectral_axis is not None:
            # Ensure that the spectral axis is an astropy Quantity
            if not isinstance(spectral_axis, u.Quantity):
                raise ValueError("Spectral axis must be a `Quantity` or "
                                 "`SpectralAxis` object.")

            # If spectral axis is provided as an astropy Quantity, convert it
            # to a specutils SpectralAxis object.
            if not isinstance(spectral_axis, SpectralAxis):
                if spectral_axis.shape[0] == flux.shape[-1] + 1:
                    bin_specification = "edges"
                else:
                    bin_specification = "centers"
                self._spectral_axis = SpectralAxis(
                    spectral_axis, redshift=redshift,
                    radial_velocity=radial_velocity, doppler_rest=rest_value,
                    doppler_convention=velocity_convention,
                    bin_specification=bin_specification)
            # If a SpectralAxis object is provided, we assume it doesn't need
            # information from other keywords added
            else:
                for a in [radial_velocity, redshift]:
                    if a is not None:
                        raise ValueError("Cannot separately set redshift or "
                                         "radial_velocity if a SpectralAxis "
                                         "object is input to spectral_axis")

                self._spectral_axis = spectral_axis

            if wcs is None:
                wcs = gwcs_from_array(self._spectral_axis)

        elif wcs is None:
            # If no spectral axis or wcs information is provided, initialize
            # with an empty gwcs based on the flux.
            size = len(flux) if not flux.isscalar else 1
            wcs = gwcs_from_array(np.arange(size) * u.Unit(""))

        super().__init__(
            data=flux.value if isinstance(flux, u.Quantity) else flux,
            wcs=wcs, **kwargs
            )

        # If no spectral_axis was provided, create a SpectralCoord based on
        # the WCS
        if spectral_axis is None:
            # If the WCS doesn't have a spectral attribute, we assume it's the
            # dummy GWCS we created or a solely spectral WCS
            if hasattr(self.wcs, "spectral"):
                # Handle generated 1D WCS that aren't set to spectral
                if not self.wcs.is_spectral and self.wcs.naxis == 1:
                    spec_axis = self.wcs.pixel_to_world(np.arange(self.flux.shape[-1]))
                else:
                    spec_axis = self.wcs.spectral.pixel_to_world(np.arange(self.flux.shape[-1]))
            else:
                spec_axis = self.wcs.pixel_to_world(np.arange(self.flux.shape[-1]))

            try:
                if spec_axis.unit.is_equivalent(u.one):
                    spec_axis = spec_axis * u.pixel
            except AttributeError:
                raise AttributeError(f"spec_axis does not have unit: "
                                     f"{type(spec_axis)} {spec_axis}")

            self._spectral_axis = SpectralAxis(
                spec_axis,
                redshift=redshift, radial_velocity=radial_velocity,
                doppler_rest=rest_value,
                doppler_convention=velocity_convention)

        if hasattr(self, 'uncertainty') and self.uncertainty is not None:
            if not flux.shape == self.uncertainty.array.shape:
                raise ValueError(
                    "Flux axis ({}) and uncertainty ({}) shapes must be the "
                    "same.".format(flux.shape, self.uncertainty.array.shape))

    def __getitem__(self, item):
        """
        Override the class indexer. We do this here because there are two cases
        for slicing on a ``Spectrum1D``:

            1.) When the flux is one dimensional, indexing represents a single
                flux value at a particular spectral axis bin, and returns a new
                ``Spectrum1D`` where all attributes are sliced.
            2.) When flux is multi-dimensional (i.e. several fluxes over the
                same spectral axis), indexing returns a new ``Spectrum1D`` with
                the sliced flux range and a deep copy of all other attributes.

        The first case is handled by the parent class, while the second is
        handled here.
        """

        if self.flux.ndim > 1 or (type(item) == tuple and item[0] == Ellipsis):
            if type(item) == tuple:
                if len(item) == len(self.flux.shape) or item[0] == Ellipsis:
                    spec_item = item[-1]
                    if not isinstance(spec_item, slice):
                        if isinstance(item, u.Quantity):
                            raise ValueError("Indexing on single spectral axis "
                                             "values is not currently allowed, "
                                             "please use a slice.")
                        spec_item = slice(spec_item, spec_item+1, None)
                        item = item[:-1] + (spec_item,)
                else:
                    # Slicing on less than the full number of axes means we want
                    # to keep the whole spectral axis
                    spec_item = slice(None, None, None)
            elif isinstance(item, slice) and (isinstance(item.start, u.Quantity) or
                    isinstance(item.stop, u.Quantity)):
                # We only allow slicing with world coordinates along the spectral
                # axis for now
                for attr in ("start", "stop"):
                    if getattr(item, attr) is None:
                        continue
                    if not getattr(item, attr).unit.is_equivalent(u.AA,
                            equivalencies=u.spectral()):
                        raise ValueError("Slicing with world coordinates is only"
                                         " enabled for spectral coordinates.")
                        break
                spec_item = item
            else:
                # Slicing with a single integer or slice uses the leading axis,
                # so we keep the whole spectral axis, which is last
                spec_item = slice(None, None, None)

            if (isinstance(spec_item.start, u.Quantity) or
                    isinstance(spec_item.stop, u.Quantity)):
                temp_spec = self._spectral_slice(spec_item)
                if spec_item is item:
                    return temp_spec
                else:
                    # Drop the spectral axis slice and perform only the spatial part
                    return temp_spec[item[:-1]]

            return self._copy(
                flux=self.flux[item],
                spectral_axis=self.spectral_axis[spec_item],
                uncertainty=self.uncertainty[item]
                if self.uncertainty is not None else None,
                mask=self.mask[item] if self.mask is not None else None)

        if not isinstance(item, slice):
            if isinstance(item, u.Quantity):
                raise ValueError("Indexing on a single spectral axis values is not"
                                 " currently allowed, please use a slice.")
            # Handle tuple slice as input by NDCube crop method
            elif isinstance(item, tuple):
                if len(item) == 1 and isinstance(item[0], slice):
                    item = item[0]
                else:
                    raise ValueError(f"Unclear how to slice with tuple {item}")
            else:
                item = slice(item, item + 1, None)
        elif (isinstance(item.start, u.Quantity) or isinstance(item.stop, u.Quantity)):
            return self._spectral_slice(item)

        tmp_spec = super().__getitem__(item)

        # TODO: this is a workaround until we figure out how to deal with non-
        #  strictly ascending spectral axes. Currently, the wcs is created from
        #  a spectral axis array by converting to a length physical type. On
        #  a regular slicing operation, the wcs is handed back to the
        #  initializer and a new spectral axis is created. This would then also
        #  be in length units, which may not be the units used initially. So,
        #  we create a new ``Spectrum1D`` that includes the sliced spectral
        #  axis. This means that a new wcs object will be created with the
        #  appropriate unit translation handling.
        return tmp_spec._copy(
            spectral_axis=self.spectral_axis[item])

    def _copy(self, **kwargs):
        """
        Perform deep copy operations on each attribute of the ``Spectrum1D``
        object.
        """
        alt_kwargs = dict(
            flux=deepcopy(self.flux),
            spectral_axis=deepcopy(self.spectral_axis),
            uncertainty=deepcopy(self.uncertainty),
            wcs=deepcopy(self.wcs),
            mask=deepcopy(self.mask),
            meta=deepcopy(self.meta),
            unit=deepcopy(self.unit),
            velocity_convention=deepcopy(self.velocity_convention),
            rest_value=deepcopy(self.rest_value))

        alt_kwargs.update(kwargs)

        return self.__class__(**alt_kwargs)

    def _spectral_slice(self, item):
        """
        Perform a region extraction given a slice on the spectral axis.
        """
        from ..manipulation import extract_region
        if item.start is None:
            start = self.spectral_axis[0]
        else:
            start = item.start
        if item.stop is None:
            stop = self.spectral_axis[-1]
        else:
            # Force the upper bound to be open, as in normal python array slicing
            exact_match = np.where(self.spectral_axis == item.stop)
            if len(exact_match[0]) == 1:
                stop_index = exact_match[0][0] - 1
                stop = self.spectral_axis[stop_index]
            else:
                stop = item.stop
        reg = SpectralRegion(start, stop)
        return extract_region(self, reg)

    @NDCube.mask.setter
    def mask(self, value):
        # Impose stricter checks than the base NDData mask setter
        if value is not None:
            value = np.array(value)
            if not self.data.shape == value.shape:
                raise ValueError(
                    "Flux axis ({}) and mask ({}) shapes must be the "
                    "same.".format(self.data.shape, value.shape))
        self._mask = value

    @property
    def frequency(self):
        """
        The `spectral_axis` as a `~astropy.units.Quantity` in units of GHz
        """
        return self.spectral_axis.to(u.GHz, u.spectral())

    @property
    def wavelength(self):
        """
        The `spectral_axis` as a `~astropy.units.Quantity` in units of Angstroms
        """
        return self.spectral_axis.to(u.AA, u.spectral())

    @property
    def energy(self):
        """
        The energy of the spectral axis as a `~astropy.units.Quantity` in units
        of eV.
        """
        return self.spectral_axis.to(u.eV, u.spectral())

    @property
    def photon_flux(self):
        """
        The flux density of photons as a `~astropy.units.Quantity`, in units of
        photons per cm^2 per second per spectral_axis unit
        """
        flux_in_spectral_axis_units = self.flux.to(
            u.W * u.cm**-2 * self.spectral_axis.unit**-1,
            u.spectral_density(self.spectral_axis))
        photon_flux_density = flux_in_spectral_axis_units / (self.energy / u.photon)

        return photon_flux_density.to(u.photon * u.cm**-2 * u.s**-1 *
                                      self.spectral_axis.unit**-1)

    @lazyproperty
    def bin_edges(self):
        return self.spectral_axis.bin_edges

    @property
    def shape(self):
        return self.flux.shape

    @property
    def redshift(self):
        """
        The redshift(s) of the objects represented by this spectrum.  May be
        scalar (if this spectrum's ``flux`` is 1D) or vector.  Note that
        the concept of "redshift of a spectrum" can be ambiguous, so the
        interpretation is set to some extent by either the user, or operations
        (like template fitting) that set this attribute when they are run on
        a spectrum.
        """
        return self.spectral_axis.redshift

    @redshift.setter
    def redshift(self, val):
        new_spec_coord = self.spectral_axis.with_radial_velocity_shift(
            -self.spectral_axis.radial_velocity).with_radial_velocity_shift(val)
        self._spectral_axis = new_spec_coord

    @property
    def radial_velocity(self):
        """
        The radial velocity(s) of the objects represented by this spectrum.  May
        be scalar (if this spectrum's ``flux`` is 1D) or vector.  Note that
        the concept of "RV of a spectrum" can be ambiguous, so the
        interpretation is set to some extent by either the user, or operations
        (like template fitting) that set this attribute when they are run on
        a spectrum.
        """
        return self.spectral_axis.radial_velocity

    @radial_velocity.setter
    def radial_velocity(self, val):
        if val is not None:
            if not val.unit.is_equivalent(u.km/u.s):
                raise u.UnitsError("Radial velocity must be a velocity.")

        new_spectral_axis = self.spectral_axis.with_radial_velocity_shift(
            -self.spectral_axis.radial_velocity).with_radial_velocity_shift(val)
        self._spectral_axis = new_spectral_axis

    def __add__(self, other):
        if not isinstance(other, NDCube):
            other = u.Quantity(other, unit=self.unit)

        return self.add(other)

    def __sub__(self, other):
        if not isinstance(other, NDCube):
            other = u.Quantity(other, unit=self.unit)

        return self.subtract(other)

    def __mul__(self, other):
        if not isinstance(other, NDCube):
            other = u.Quantity(other)

        return self.multiply(other)

    def __div__(self, other):
        if not isinstance(other, NDCube):
            other = u.Quantity(other)

        return self.divide(other)

    def __truediv__(self, other):
        if not isinstance(other, NDCube):
            other = u.Quantity(other)

        return self.divide(other)

    def _format_array_summary(self, label, array):
        if len(array) == 1:
            mean = np.mean(array)
            s = "{:17} [ {:.5} ],  mean={:.5}"
            return s.format(label+':', array[0], array[-1], mean)
        elif len(array) > 1:
            mean = np.mean(array)
            s = "{:17} [ {:.5}, ..., {:.5} ],  mean={:.5}"
            return s.format(label+':', array[0], array[-1], mean)
        else:
            return "{:17} [ ],  mean= n/a".format(label+':')

    def __str__(self):
        result = "Spectrum1D "
        # Handle case of single value flux
        if self.flux.ndim == 0:
            result += "(length=1)\n"
            return result + "flux:   {}".format(self.flux)

        # Handle case of multiple flux arrays
        result += "(length={})\n".format(len(self.spectral_axis))
        if self.flux.ndim > 1:
            for i, flux in enumerate(self.flux):
                label = 'flux{:2}'.format(i)
                result += self._format_array_summary(label, flux) + '\n'
        else:
            result += self._format_array_summary('flux', self.flux) + '\n'
        # Add information about spectral axis
        result += self._format_array_summary('spectral axis', self.spectral_axis)
        # Add information about uncertainties if available
        if self.uncertainty:
            result += "\nuncertainty:      [ {}, ..., {} ]".format(
                self.uncertainty[0], self.uncertainty[-1])
        return result

    def __repr__(self):
        inner_str = "flux={}, spectral_axis={}".format(repr(self.flux),
                                                       repr(self.spectral_axis))

        if self.uncertainty is not None:
            inner_str += ", uncertainty={}".format(repr(self.uncertainty))

        result = "<Spectrum1D({})>".format(inner_str)

        return result
