import warnings
from copy import deepcopy

import numpy as np
from astropy import units as u
from astropy.coordinates import SpectralCoord
from astropy.utils.decorators import lazyproperty
from astropy.utils.decorators import deprecated
from astropy.nddata import NDUncertainty, NDIOMixin, NDArithmeticMixin, NDDataArray
from gwcs.wcs import WCS as GWCS

from .spectral_axis import SpectralAxis
from .spectrum_mixin import OneDSpectrumMixin
from .spectral_region import SpectralRegion
from ..utils.wcs_utils import gwcs_from_array

from ndcube import NDCube

__all__ = ['Spectrum1D', 'Spectrum']


class Spectrum(OneDSpectrumMixin, NDCube, NDIOMixin, NDArithmeticMixin):
    """
    Spectrum container for N-dimensional data with one spectral axis.

    `Spectrum` can contain "vector spectra" by having the
    ``flux`` have a shape with dimension greater than 1.  One dimension of
    ``flux`` must match the length of the ``spectral_axis`` if ``spectral_axis``
    is provided.

    For multidimensional spectra that are all the same shape but have different
    spectral axes, use a :class:`~specutils.SpectrumCollection`.  For a
    collection of spectra that have different shapes, use
    :class:`~specutils.SpectrumList`. For more on this topic, see
    :ref:`specutils-representation-overview`.

    Parameters
    ----------
    flux : `~astropy.units.Quantity`
        The flux data for this spectrum. This can be a simple `~astropy.units.Quantity`,
        or an existing `~Spectrum` or `~ndcube.NDCube` object.
    spectral_axis : `~astropy.units.Quantity` or `~specutils.SpectralAxis`
        Dispersion information with the same shape as the last (or only)
        dimension of flux, or one greater than the last dimension of flux
        if specifying bin edges.
    spectral_axis_index : integer, optional
        If it is ambiguous which axis is the spectral axis (e.g., if there are multiple
        axes in the flux array with the same length as the input spectral_axis),
        this argument is used to specify which is the spectral axis.
    move_spectral_axis : int, str, optional
        Force the spectral axis to be either the last axis (the default behavior prior
        to version 2.0) by setting this argument to 'last' or -1, or the first axis by
        setting this argument to 'first' or 0.
        This will do a simple ``swapaxis`` between the relevant axis and original
        spectral axis. If None, the spectral axis is left wherever it is in the input.
    wcs : `~astropy.wcs.WCS` or `~gwcs.wcs.WCS`
        WCS information object that either has a spectral component or is
        only spectral.
    velocity_convention : {"relativistic", "optical", "radio"}
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
    def __init__(self, flux=None, spectral_axis=None, spectral_axis_index=None,
                 wcs=None, velocity_convention=None, rest_value=None,
                 redshift=None, radial_velocity=None, bin_specification=None,
                 move_spectral_axis=None, **kwargs):

        if spectral_axis_index == -1:
            spectral_axis_index = flux.ndim - 1

        # If the flux (data) argument is already a Spectrum (as it would
        # be for internal arithmetic operations), avoid setup entirely.
        if isinstance(flux, Spectrum):
            self._spectral_axis_index = flux.spectral_axis_index
            self._spectral_axis = flux.spectral_axis
            super().__init__(flux)
            return

        self._spectral_axis_index = spectral_axis_index
        # Might as well handle this right away
        if spectral_axis_index is None and flux is not None:
            if flux.ndim == 1:
                self._spectral_axis_index = 0
        elif flux is None:
            self._spectral_axis_index = 0

        # Check for pre-defined entries in the kwargs dictionary.
        unknown_kwargs = set(kwargs).difference(
            {'data', 'unit', 'uncertainty', 'meta', 'mask', 'copy',
             'extra_coords'})

        if len(unknown_kwargs) > 0:
            raise ValueError("Initializer contains unknown arguments(s): {}."
                             "".format(', '.join(map(str, unknown_kwargs))))

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
        if np.ndim(flux) == 0 and spectral_axis is None and wcs is None:
            super(Spectrum, self).__init__(data=flux, wcs=wcs, **kwargs)
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
                warnings.warn("No unit information provided with rest value. "
                              f"Assuming units of spectral axis ('{spectral_axis.unit}').")
                rest_value = u.Quantity(rest_value, spectral_axis.unit)
            elif not rest_value.unit.is_equivalent(u.AA, equivalencies=u.spectral()):
                raise u.UnitsError("Rest value must be "
                                   "energy/wavelength/frequency equivalent.")

        # If flux and spectral axis are both specified, check that their lengths
        # match or are off by one (implying the spectral axis stores bin edges).
        # If we can't determine which flux axis corresponds to the spectral axis
        # we raise an error.
        if flux is not None and spectral_axis is not None:
            if spectral_axis_index is None:
                if flux.ndim == 1:
                    self._spectral_axis_index = 0
                else:
                    matching_axes = []
                    if bin_specification == "centers":
                        add_elements = [0,]
                    elif bin_specification == "edges":
                        add_elements = [1,]
                    elif bin_specification is None:
                        add_elements = [0,1]
                    for i in range(flux.ndim):
                        for add_element in add_elements:
                            if spectral_axis.shape[0] == flux.shape[i] + add_element:
                                matching_axes.append(i)

                    if len(matching_axes) == 1:
                        self._spectral_axis_index = matching_axes[0]
                    else:
                        raise ValueError("Unable to determine which flux axis corresponds to "
                                         "the spectral axis. Please specify spectral_axis_index"
                                         " or provide a spectral_axis matching a flux axis.")

            # Make sure the length of the spectral axis matches the appropriate flux axis
            if spectral_axis.shape[0] == flux.shape[self.spectral_axis_index]:
                if bin_specification == "edges":
                    raise ValueError("A spectral axis input as bin edges must "
                                     "have length one greater than the flux axis")
                bin_specification = "centers"
            elif spectral_axis.shape[0] == flux.shape[self.spectral_axis_index]+1:
                if bin_specification == "centers":
                    raise ValueError("A spectral axis input as bin centers "
                        "must be the same length as the flux axis")
                bin_specification = "edges"
            else:
                raise ValueError(
                    f"Spectral axis length ({spectral_axis.shape[0]}) must be the "
                    "same size or one greater (if specifying bin edges) than that "
                    f"of the corresponding flux axis ({flux.shape[self.spectral_axis_index]})")

        # If a WCS is provided, determine which axis is the spectral axis
        if wcs is not None:
            naxis = None
            if hasattr(wcs, "naxis"):
                naxis = wcs.naxis
            # GWCS doesn't have naxis
            elif hasattr(wcs, "world_n_dim"):
                naxis = wcs.world_n_dim

            if naxis is not None and naxis > 1:
                temp_axes = []
                phys_axes = wcs.world_axis_physical_types
                if self._spectral_axis_index is None:
                    for i in range(len(phys_axes)):
                        if phys_axes[i] is None:
                            continue
                        if (phys_axes[i][0:2] == "em" or phys_axes[i][0:5] == "spect" or
                                phys_axes[i][7:12] == "Spect"):
                            temp_axes.append(i)
                    if len(temp_axes) != 1:
                        raise ValueError("Input WCS must have exactly one axis with "
                                        "spectral units, found {}".format(len(temp_axes)))
                    else:
                        # Due to FITS conventions, the WCS axes are listed in opposite
                        # order compared to the data array.
                        self._spectral_axis_index = len(flux.shape)-temp_axes[0]-1

                if move_spectral_axis is not None:
                    if isinstance(wcs, GWCS):
                        raise ValueError("move_spectral_axis cannot be used with GWCS")
                    if isinstance(move_spectral_axis, str):
                        if move_spectral_axis.lower() == 'first':
                            move_to_index = 0
                        elif move_spectral_axis.lower() == 'last':
                            move_to_index = wcs.naxis - 1
                        else:
                            raise ValueError("move_spectral_axis must be either 'first' or 'last'")
                    elif isinstance(move_spectral_axis, int):
                        move_to_index = move_spectral_axis
                    else:
                        raise ValueError("move_spectral_axis must be an integer or 'first'/'last'")

                    if move_to_index != self._spectral_axis_index:
                        wcs = wcs.swapaxes(self._spectral_axis_index, move_to_index)
                        if flux is not None:
                            flux = np.swapaxes(flux, self._spectral_axis_index, move_to_index)
                        if "mask" in kwargs:
                            if kwargs["mask"] is not None:
                                kwargs["mask"] = np.swapaxes(kwargs["mask"],
                                                    self._spectral_axis_index, move_to_index)
                        if "uncertainty" in kwargs:
                            if kwargs["uncertainty"] is not None:
                                if isinstance(kwargs["uncertainty"], NDUncertainty):
                                    # Account for Astropy uncertainty types
                                    temp_unc = np.swapaxes(kwargs["uncertainty"].array,
                                                           self._spectral_axis_index, move_to_index)
                                    if kwargs["uncertainty"].unit is not None:
                                        temp_unc = temp_unc * u.Unit(kwargs["uncertainty"].unit)
                                    kwargs["uncertainty"] = type(kwargs["uncertainty"])(temp_unc)
                                else:
                                    kwargs["uncertainty"] = np.swapaxes(kwargs["uncertainty"],
                                                            self._spectral_axis_index, move_to_index)

                        self._spectral_axis_index = move_to_index

            else:
                if flux is not None and flux.ndim == 1:
                    self._spectral_axis_index = 0
                else:
                    if self.spectral_axis_index is None:
                        raise ValueError("WCS is 1D but flux is multi-dimensional. Please"
                                         " specify spectral_axis_index.")

        elif move_spectral_axis is not None:
            raise ValueError("Unable to use `move_spectral_axis` without a multi-dimensional WCS")

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
                wcs = gwcs_from_array(self._spectral_axis,
                                      flux.shape,
                                      spectral_axis_index=self.spectral_axis_index
                                      )

        elif wcs is None:
            # If no spectral axis or wcs information is provided, initialize
            # with an empty gwcs based on the flux.
            if self.spectral_axis_index is None:
                if flux.ndim == 1:
                    self._spectral_axis_index = 0
                else:
                    raise ValueError("Must specify spectral_axis_index if no WCS or spectral"
                                     " axis is input.")
            size = flux.shape[self.spectral_axis_index] if not flux.isscalar else 1
            wcs = gwcs_from_array(np.arange(size) * u.Unit("pixel"),
                                  flux.shape,
                                  spectral_axis_index=self.spectral_axis_index
                                  )

        super().__init__(
            data=flux.value if isinstance(flux, u.Quantity) else flux,
            wcs=wcs, **kwargs)

        # If no spectral_axis was provided, create a SpectralCoord based on
        # the WCS
        if spectral_axis is None:
            # If the WCS doesn't have a spectral attribute, we assume it's the
            # dummy GWCS we created or a solely spectral WCS
            if hasattr(self.wcs, "spectral"):
                # Handle generated 1D WCS that aren't set to spectral
                if not self.wcs.is_spectral and self.wcs.naxis == 1:
                    spec_axis = self.wcs.pixel_to_world(
                                    np.arange(self.flux.shape[self.spectral_axis_index]))
                else:
                    spec_axis = self.wcs.spectral.pixel_to_world(
                                    np.arange(self.flux.shape[self.spectral_axis_index]))
            else:
                # We now keep the entire GWCS, including spatial information, so we need to include
                # all axes in the pixel_to_world call. Note that this assumes/requires that the
                # dispersion is the same at all spatial locations.
                wcs_args = []
                for i in range(len(self.flux.shape)):
                    wcs_args.append(np.zeros(self.flux.shape[self.spectral_axis_index]))
                # Replace with arange for the spectral axis
                wcs_args[self.spectral_axis_index] = np.arange(self.flux.shape[self.spectral_axis_index])
                wcs_args.reverse()
                temp_coords = self.wcs.pixel_to_world(*wcs_args)
                # If there are spatial axes, temp_coords will have a SkyCoord and a SpectralCoord
                if isinstance(temp_coords, list):
                    for coords in temp_coords:
                        if isinstance(coords, SpectralCoord):
                            spec_axis = coords
                            break
                    else:
                        # WCS axis ordering is reverse of numpy
                        spec_axis = temp_coords[len(temp_coords) - self.spectral_axis_index - 1]
                else:
                    spec_axis = temp_coords

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

        # make sure that spectral axis is strictly increasing or decreasing,
        # raise an error if not.
        if not self._check_strictly_increasing_decreasing():
            raise ValueError('Spectral axis must be strictly increasing or decreasing.')

        if hasattr(self, 'uncertainty') and self.uncertainty is not None:
            if not flux.shape == self.uncertainty.array.shape:
                raise ValueError(
                    "Flux axis ({}) and uncertainty ({}) shapes must be the "
                    "same.".format(flux.shape, self.uncertainty.array.shape))

    def __getitem__(self, item):
        """
        Override the class indexer. We do this here because there are two cases
        for slicing on a ``Spectrum``:

            1.) When the flux is one dimensional, indexing represents a single
                flux value at a particular spectral axis bin, and returns a new
                ``Spectrum`` where all attributes are sliced.
            2.) When flux is multi-dimensional (i.e. several fluxes over the
                same spectral axis), indexing returns a new ``Spectrum`` with
                the sliced flux range and a deep copy of all other attributes.

        The first case is handled by the parent class, while the second is
        handled here.
        """
        new_spectral_axis_index = self.spectral_axis_index

        if self.flux.ndim > 1 or (isinstance(item, tuple) and item[0] is Ellipsis):
            if isinstance(item, tuple):
                if len(item) == self.flux.ndim or item[0] == Ellipsis:
                    spec_item = item[self.spectral_axis_index]
                    if not isinstance(spec_item, slice):
                        if isinstance(item, u.Quantity):
                            raise ValueError("Indexing on single spectral axis "
                                             "values is not currently allowed, "
                                             "please use a slice.")
                        spec_item = slice(spec_item, spec_item+1, None)
                        # We have to hack around updating this tuple entry
                        item = list(item)
                        item[self.spectral_axis_index] = spec_item
                        item = tuple(item)
                else:
                    # Slicing on less than the full number of axes means we want
                    # to keep the whole spectral axis
                    spec_item = slice(None, None, None)
                    # If any slices are single integers, need to decrement the spectral axis index

                for i in range(len(item)-1):
                    # Decrement spectral_axis_index for each single element slice
                    if isinstance(item[i], int):
                        new_spectral_axis_index -= 1

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
                # so we keep the whole spectral axis, which is last. Need to decrement
                # the spectral axis index by one
                if isinstance(item, int):
                    new_spectral_axis_index -= 1
                spec_item = slice(None, None, None)

            if (isinstance(spec_item.start, u.Quantity) or
                    isinstance(spec_item.stop, u.Quantity)):
                temp_spec = self._spectral_slice(spec_item)
                if spec_item is item:
                    return temp_spec
                else:
                    # Drop the spectral axis slice and perform only the spatial part
                    return temp_spec[item[:-1]]

            if "original_wcs" not in self.meta:
                new_meta = deepcopy(self.meta)
                new_meta["original_wcs"] = deepcopy(self.wcs)
            else:
                new_meta = deepcopy(self.meta)

            return self._copy(
                flux=self.flux[item],
                spectral_axis=self.spectral_axis[spec_item],
                uncertainty=self.uncertainty[item]
                if self.uncertainty is not None else None,
                mask=self.mask[item] if self.mask is not None else None,
                meta=new_meta, wcs=None, spectral_axis_index=new_spectral_axis_index)

        if not isinstance(item, slice):
            if isinstance(item, u.Quantity):
                raise ValueError("Indexing on a single spectral axis value is not"
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

        # Work around error in SpectralCoord creation in super().__getitem__ for
        # spectral axis in pixels.
        if self.spectral_axis.unit == u.pix:
            if "original_wcs" not in self.meta:
                new_meta = deepcopy(self.meta)
                new_meta["original_wcs"] = deepcopy(self.wcs)
            return self._copy(
                flux=self.flux[item],
                spectral_axis=self.spectral_axis[item],
                uncertainty=self.uncertainty[item]
                if self.uncertainty is not None else None,
                mask=self.mask[item] if self.mask is not None else None,
                meta=new_meta, wcs=None, spectral_axis_index=new_spectral_axis_index)

        tmp_spec = super().__getitem__(item)

        # TODO: this is a workaround until we figure out how to deal with non-
        #  strictly ascending spectral axes. Currently, the wcs is created from
        #  a spectral axis array by converting to a length physical type. On
        #  a regular slicing operation, the wcs is handed back to the
        #  initializer and a new spectral axis is created. This would then also
        #  be in length units, which may not be the units used initially. So,
        #  we create a new ``Spectrum`` that includes the sliced spectral
        #  axis. This means that a new wcs object will be created with the
        #  appropriate unit translation handling.
        if "original_wcs" not in self.meta:
            new_meta = deepcopy(self.meta)
            new_meta["original_wcs"] = deepcopy(self.wcs)
        else:
            new_meta = deepcopy(self.meta)

        return tmp_spec._copy(spectral_axis=self.spectral_axis[item], wcs=None,
                              meta=new_meta)

    def _copy(self, **kwargs):
        """
        Perform deep copy operations on each attribute of the ``Spectrum``
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
            rest_value=deepcopy(self.rest_value),
            spectral_axis_index=deepcopy(self.spectral_axis_index))

        alt_kwargs.update(kwargs)
        if 'spectral_axis' in kwargs and 'wcs' not in kwargs:
            # We assume in this case that the user wants to override with a new spectral axis
            alt_kwargs['meta']['original_wcs'] = alt_kwargs['wcs']
            alt_kwargs['wcs'] = None

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

    def collapse(self, method, axis=None):
        """
        Collapse the flux array given a method. Will collapse either to a single
        value (default), over a specified numerical axis or axes if specified, or
        over the spectral or non-spectral axes if ``physical_type`` is specified.

        If the collapse leaves the spectral axis unchanged, a `~specutils.Spectrum`
        will be returned. Otherwise an `~astropy.units.Quantity` array will be
        returned.

        Note that these calculations are not currently uncertainty-aware, but do
        respect masks.

        Parameters
        ----------

        method : str, function
            The method by which the flux will be collapsed. String options are
            'mean', 'min', 'max', 'sum', and 'median'. Also accepts a function
            as input, which must take an `astropy.units.Quantity` array as input
            and accept an 'axis' argument.
        axis : int, tuple, str, optional
            The axis or axes over which to collapse the flux array. May also be
            a string, either 'spectral' to collapse over the spectral axis, or
            'spatial' to collapse over all other axes.

        Returns
        -------
        :class:`~specutils.Spectrum` or :class:`~astropy.units.Quantity`

        """
        collapse_funcs = {"mean": np.nanmean, "max": np.nanmax, "min": np.nanmin,
                         "median": np.nanmedian, "sum": np.nansum}

        if isinstance(axis, str):
            if axis == 'spectral':
                axis = self.spectral_axis_index
            elif axis == 'spatial':
                # generate tuple if needed for multiple spatial axes
                axis = tuple([x for x in range(len(self.flux.shape)) if
                              x != self.spectral_axis_index])
            else:
                raise ValueError("String axis input must be 'spatial' or 'spectral'")

        # Set masked locations to NaN for the calculation, since the `where` argument
        # does not seem to work consistently in the numpy functions.
        flux_to_collapse = self.flux.copy()
        if self.mask is not None:
            flux_to_collapse[np.where(self.mask != 0)] = np.nan

        # Leave open the possibility of the user providing their own method
        if callable(method):
            collapsed_flux = method(flux_to_collapse, axis=axis)
        else:
            collapsed_flux = collapse_funcs[method](flux_to_collapse, axis=axis)

        # Return a Spectrum if we collapsed over the spectral axis, a Quantity if not
        if axis in (self.spectral_axis_index, None):
            return collapsed_flux
        elif isinstance(axis, tuple) and self.spectral_axis_index in axis:
            return collapsed_flux
        else:
            # Pass the spectral axis rather than WCS in this case, so we don't have to
            # figure out which part of a multidimensional WCS is the spectral part.
            return Spectrum(collapsed_flux, spectral_axis=self.spectral_axis)

    def mean(self, **kwargs):
        return self.collapse("mean", **kwargs)

    def max(self, **kwargs):
        return self.collapse("max", **kwargs)

    def min(self, **kwargs):
        return self.collapse("min", **kwargs)

    def median(self, **kwargs):
        return self.collapse("median", **kwargs)

    def sum(self, **kwargs):
        return self.collapse("sum", **kwargs)

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
    def spectral_axis_direction(self):
        return self._spectral_axis_direction

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

    def set_redshift_to(self, redshift):
        """
        This sets the redshift of the spectrum to be `redshift` *without*
        changing the values of the `spectral_axis`.

        If you want to shift the `spectral_axis` based on this value, use
        `shift_spectrum_to`.
        """
        new_spec_coord = self.spectral_axis.replicate(redshift=redshift)
        self._spectral_axis = new_spec_coord

    def set_radial_velocity_to(self, radial_velocity):
        """
        This sets the radial velocity of the spectrum to be `radial_velocity`
        *without* changing the values of the `spectral_axis`.

        If you want to shift the `spectral_axis` based on this value, use
        `shift_spectrum_to`.
        """
        new_spec_coord = self.spectral_axis.replicate(
            radial_velocity=radial_velocity
        )
        self._spectral_axis = new_spec_coord

    def shift_spectrum_to(self, *, redshift=None, radial_velocity=None):
        """
        This shifts in-place the values of the `spectral_axis`, given either a
        redshift or radial velocity.

        If you do *not* want to change the `spectral_axis`, use
        `set_redshift_to` or `set_radial_velocity_to`.
        """
        if redshift is not None and radial_velocity is not None:
            raise ValueError(
                "Only one of redshift or radial_velocity can be used."
            )
        if redshift is not None:
            new_spec_coord = self.spectral_axis.with_radial_velocity_shift(
                -self.spectral_axis.radial_velocity
            ).with_radial_velocity_shift(redshift)
            self._spectral_axis = new_spec_coord
        elif radial_velocity is not None:
            if radial_velocity is not None:
                if not radial_velocity.unit.is_equivalent(u.km/u.s):
                    raise u.UnitsError("Radial velocity must be a velocity.")

            new_spectral_axis = self.spectral_axis.with_radial_velocity_shift(
                -self.spectral_axis.radial_velocity
            ).with_radial_velocity_shift(radial_velocity)
            self._spectral_axis = new_spectral_axis
        else:
            raise ValueError("One of redshift or radial_velocity must be set.")

    def with_spectral_axis_last(self):
        """
        Convenience method to return a new copy of the Spectrum with the spectral axis last.
        """
        return Spectrum(flux=self.flux, wcs=self.wcs,
                          mask=self.mask, uncertainty=self.uncertainty,
                          redshift=self.redshift, move_spectral_axis="last")

    def _return_with_redshift(self, result):
        # We need actual spectral units to shift
        if result.spectral_axis.unit not in ('', 'pix', 'pixels'):
            result.shift_spectrum_to(redshift=self.redshift)
        return result

    def _check_input(self, other, force_quantity=False):
        # NDArithmetic mixin will try to turn other into a Spectrum, which will fail
        # sometimes because of not specifiying the spectral axis index
        if isinstance(other, Spectrum):
            # Take this opportunity to check if the spectral axes match
            if not np.all(other.spectral_axis == self.spectral_axis):
                raise ValueError("Spectral axis of both operands must match")
        else:
            if not isinstance(other, u.Quantity) and force_quantity:
                other = other * self.unit

        return other

    def _do_flux_arithmetic(self, other, arith_func):
        '''
        Perform an arithmetic operation by casting the flux as a NDDataArray
        '''
        operand1 = NDDataArray(self.flux, uncertainty=self.uncertainty, mask=self.mask)
        if isinstance(other, (Spectrum)):
            other = NDDataArray(other.flux, uncertainty=other.uncertainty, mask=other.mask)

        func = getattr(operand1, arith_func)
        new_flux = func(other)
        return self._return_with_redshift(Spectrum(new_flux.data*new_flux.unit,
                                                   wcs=self.wcs,
                                                   meta=self.meta,
                                                   uncertainty=new_flux.uncertainty,
                                                   mask = new_flux.mask,
                                                   spectral_axis_index=self.spectral_axis_index))

    def __add__(self, other):
        other = self._check_input(other, force_quantity=True)
        return self._do_flux_arithmetic(other, "add")

    def __sub__(self, other):
        try:
            other = self._check_input(other, force_quantity=True)
        except TypeError:
            # Might need special handling by other operand before ours
            if hasattr(other, "__rsub__"):
                return other.__rsub__(self)
            else:
                raise

        return self._do_flux_arithmetic(other, "subtract")

    def __mul__(self, other):
        other = self._check_input(other)
        return self._do_flux_arithmetic(other, "multiply")

    def __div__(self, other):
        other = self._check_input(other)
        return self._do_flux_arithmetic(other, "divide")

    def __truediv__(self, other):
        other = self._check_input(other)
        return self._do_flux_arithmetic(other, "divide")

    __radd__ = __add__

    __rmul__ = __mul__

    def __rsub__(self, other):
        return -1 * (self - other)

    def _format_array_summary(self, label, array):
        array_str = np.array2string(array, threshold=8, prefix=label)
        if len(array) >= 1:
            mean = np.nanmean(array)
            s = f"{label}{array_str} {array.unit},  mean={mean:.5f}"
            return s
        else:
            return "{:17} [ ],  mean= n/a".format(label+':')

    def __str__(self):
        result = "Spectrum "
        result += "(length={})\n".format(len(self.spectral_axis))

        # Add Flux information
        result += self._format_array_summary('Flux=', self.flux) + '\n'

        # Add information about spectral axis
        result += self._format_array_summary('Spectral Axis=', self.spectral_axis)

        # Add information about uncertainties if available
        if self.uncertainty:
            result += (f'\nUncertainty={type(self.uncertainty).__name__} '
                       f'({np.array2string(self.uncertainty.array, threshold=8)}'
                       f' {self.uncertainty.unit})')

        return result

    def __repr__(self):
        flux_str  = "flux="
        if (self.flux.ndim == 1 and self.flux.size <= 10) or self.flux.size <= 20:
            flux_str += repr(self.flux)
        else:
            flux_summary = f"{self.flux.value.flat[0]} ... {self.flux.value.flat[-1]}"
            flux_str = flux_str + "[" * self.flux.ndim + flux_summary + "]" * self.flux.ndim
            flux_str += f" {self.flux.unit}"

        flux_str += f" (shape={self.flux.shape}, mean={np.nanmean(self.flux):.5f}); "
        # Sometimes this errors if an error occurs during initialization
        if hasattr(self, "_spectral_axis"):
            spectral_axis_str = (repr(self.spectral_axis).split("[")[0] +
                                    np.array2string(self.spectral_axis, threshold=8) +
                                    f" {self.spectral_axis.unit}>")
            spectral_axis_str = f"spectral_axis={spectral_axis_str} (length={len(self.spectral_axis)})"
            inner_str = (flux_str + spectral_axis_str)
        else:
            inner_str = flux_str

        if self.uncertainty is not None:
            inner_str += f"; uncertainty={self.uncertainty.__class__.__name__}"

        result = "<Spectrum({})>".format(inner_str)

        return result


@deprecated(since="2.0", alternative="Spectrum")
class Spectrum1D(Spectrum):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
