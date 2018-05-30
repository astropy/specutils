import logging

import numpy as np
from astropy import units as u
from astropy.nddata import NDDataRef
from astropy.utils.decorators import lazyproperty
from astropy.nddata import NDUncertainty
from ..wcs import WCSWrapper, WCSAdapter
from .spectrum_mixin import OneDSpectrumMixin

__all__ = ['Spectrum1D']


class Spectrum1D(OneDSpectrumMixin, NDDataRef):
    """
    Spectrum container for 1D spectral data.

    Parameters
    ----------
    flux : `numpy.ndarray`-like or `astropy.units.Quantity` or astropy.nddata.NDData`-like
        The flux data for this spectrum.
    spectral_axis : `numpy.ndarray`-like or `astropy.units.Quanitty`
        Dispersion information with the same shape as the last (or only)
        dimension of flux.
    wcs : `astropy.wcs.WCS` or `gwcs.WCS`
        WCS information object.
    unit : str or `astropy.units.Unit`
        The unit for the flux data. Must be parseable by `astropy.units.Unit`.
        If `flux` is supplied as a `astropy.units.Quantity`, this
        is superceded by the defined unit.
    spectral_axis_unit : str or `astropy.units.Unit`
        The unit for the spectral axis. Must be parseable by `astropy.units.Unit`.
        If `spectral_axis` is supplied as a `astropy.units.Quantity`, this
        is superceded by the defined unit.
    velocity_convention : {"doppler_relativistic", "doppler_optical", "doppler_radio"}
        Convention used for velocity conversions.
    rest_value : ~`astropy.units.Quantity`
        Any quantity supported by the standard spectral equivalencies
        (wavelength, energy, frequency, wave number). Describes the rest value
        of the spectral axis for use with velocity conversions.
    uncertainty : ~`astropy.nddata.NDUncertainty`
        Contains uncertainty information along with propagation rules for
        spectrum arithmetic. Can take a unit, but if none is given, will use
        the unit defined in the flux.
    meta : dict
        Arbitrary container for any user-specific information to be carried
        around with the spectrum container object.
    """

    def __init__(self, flux=None, spectral_axis=None, wcs=None, unit=None,
                 spectral_axis_unit=None, velocity_convention=None,
                 rest_value=None, *args, **kwargs):
        # In cases of slicing, new objects will be initialized with `data`
        # instead of `flux`. Ensure we grab the `data` argument.
        if flux is None and 'data' in kwargs:
            flux = kwargs.pop('data')

        # If the flux (data) argument is a subclass of nddataref (as it would
        # be for internal arithmetic operations), avoid setup entirely.
        if issubclass(flux.__class__, NDDataRef):
            self._velocity_convention = flux._velocity_convention
            self._rest_value = flux._rest_value

            super(Spectrum1D, self).__init__(flux)
            return

        # Attempt to parse the spectral axis. If none is given, try instead to
        # parse a given wcs. This is put into a GWCS object to
        # then be used behind-the-scenes for all specutils operations.
        if spectral_axis is not None:
            if not isinstance(spectral_axis, u.Quantity):
                spectral_axis = u.Quantity(spectral_axis,
                                           unit=spectral_axis_unit or u.AA)

                logging.warning("No spectral axis units given, assuming "
                                "{}".format(spectral_axis.unit))

            wcs = WCSWrapper.from_array(spectral_axis)
        elif wcs is not None:
            if not issubclass(wcs.__class__, WCSAdapter):
                wcs = WCSWrapper(wcs)
        elif isinstance(flux, float) or isinstance(flux, int) or isinstance(flux, np.ndarray):
            # In the case where the arithmetic operation is being performed with
            # a single float, int, or array object, just go ahead and ignore wcs
            # requirements
            super(Spectrum1D, self).__init__(data=flux)
            return
        else:
            # If no wcs and no spectral axis has been given, raise an error
            raise LookupError("No WCS object or spectral axis information has "
                              "been given. Please provide one.")

        if not isinstance(flux, u.Quantity):
            flux = u.Quantity(flux, unit=unit or "Jy")

        self._velocity_convention = velocity_convention

        if rest_value is None:
            if wcs.rest_frequency != 0:
                self._rest_value = wcs.rest_frequency * u.Hz
            elif wcs.rest_wavelength != 0:
                self._rest_value = wcs.rest_wavelength * u.AA
            else:
                self._rest_value = 0 * u.AA
        else:
            self._rest_value = rest_value

            if not isinstance(self._rest_value, u.Quantity):
                logging.info("No unit information provided with rest value. "
                             "Assuming units of spectral axis ('%s').",
                             spectral_axis.unit)
                self._rest_value = u.Quantity(rest_value, spectral_axis.unit)
            elif not self._rest_value.unit.is_equivalent(u.AA) and not self._rest_value.unit.is_equivalent(u.Hz):
                raise u.UnitsError("Rest value must be energy/wavelength/frequency equivalent.")

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

    @property
    def photon_flux(self):
        return (self.flux/self.energy).to(
            1.0/(self.spectral_axis.unit*u.s*(u.cm**2)))

    @lazyproperty
    def bin_edges(self):
        return self.wcs.bin_edges()

    @staticmethod
    def _compare_wcs(this_operand, other_operand):
        """
        NNData arithmetic callable to determine if two wcs's are compatible.
        """
        # If the other operand is a simple number or array, allow the operations
        if (isinstance(other_operand, float) or isinstance(other_operand, int)
            or isinstance(other_operand, np.ndarray)):
            return True

        # First check if units are equivalent, if so, create a new spectrum
        # object with spectral axis in compatible units
        other_wcs = other_operand.wcs.with_spectral_unit(
            this_operand.wcs.spectral_axis_unit)

        if other_wcs is None:
            return False

        # Check if the shape of the axes are compatible
        if this_operand.spectral_axis.shape != other_operand.spectral_axis.shape:
            logging.error("Shape of spectral axes between operands must be "
                          "equivalent.")
            return False

        # And that they cover the same range
        if (this_operand.spectral_axis[0] != other_operand.spectral_axis[0] or
                this_operand.spectral_axis[-1] != other_operand.spectral_axis[-1]):
            logging.error("Spectral axes between operands must cover the "
                          "same range. Interpolation may be required.")
            return False

        # Check if the delta dispersion is equivalent between the two axes
        if not np.array_equal(np.diff(this_operand.spectral_axis),
                              np.diff(other_operand.spectral_axis)):
            logging.error("Delta dispersion of spectral axes of operands "
                          "must be equivalent. Interpolation may be required.")
            return False

        return True

    def __add__(self, other):
        return self.add(
            other, compare_wcs=lambda o1, o2: self._compare_wcs(self, other))

    def __sub__(self, other):
        return self.subtract(
            other, compare_wcs=lambda o1, o2: self._compare_wcs(self, other))

    def __mul__(self, other):
        return self.multiply(
            other, compare_wcs=lambda o1, o2: self._compare_wcs(self, other))

    def __div__(self, other):
        return self.divide(
            other, compare_wcs=lambda o1, o2: self._compare_wcs(self, other))

    def __truediv__(self, other):
        return self.divide(
            other, compare_wcs=lambda o1, o2: self._compare_wcs(self, other))

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
