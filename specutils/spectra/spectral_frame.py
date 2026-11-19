"""
Descriptions of the medium and the velocity reference frame that a spectral
axis is expressed in.
"""
from dataclasses import dataclass, asdict

import astropy.units as u

__all__ = ['SpectralMedium', 'SPECTRAL_FRAMES', 'normalize_frame']


#: FITS ``SPECSYS`` reference frame codes (Greisen et al. 2006, Table 5),
#: mapped to the name of the astropy coordinate frame in which an observer is
#: stationary, or `None` if astropy has no equivalent frame.
SPECTRAL_FRAMES = {
    'TOPOCENT': None,  # observer moving with the telescope
    'GEOCENTR': 'gcrs',
    'BARYCENT': 'icrs',
    'HELIOCEN': 'hcrs',
    'LSRK': 'lsrk',
    'LSRD': 'lsrd',
    'GALACTOC': 'galactocentric',
    'LOCALGRP': None,
    'CMBDIPOL': None,
    'SOURCE': None,  # observer moving with the source: the rest frame
}

_FRAME_ALIASES = {
    'TOPOCENTRIC': 'TOPOCENT',
    'GEOCENTRIC': 'GEOCENTR',
    'BARYCENTRIC': 'BARYCENT',
    'HELIOCENTRIC': 'HELIOCEN',
    'GALACTOCENTRIC': 'GALACTOC',
    'LOCALGROUP': 'LOCALGRP',
    'CMBDIPOLE': 'CMBDIPOL',
    'REST': 'SOURCE',
}


def normalize_frame(frame):
    """
    Validate a spectral reference frame and return its FITS ``SPECSYS`` code.

    Parameters
    ----------
    frame : str or None
        A FITS ``SPECSYS`` code such as ``'TOPOCENT'`` or ``'BARYCENT'``
        (case-insensitive), or one of a few long-form aliases such as
        ``'barycentric'`` or ``'rest'``.

    Returns
    -------
    str or None
        The upper-case ``SPECSYS`` code, or `None` if ``frame`` is `None`.
    """
    if frame is None:
        return None
    if not isinstance(frame, str):
        raise TypeError("frame must be a string FITS SPECSYS code, "
                        f"got {type(frame).__name__}")
    key = frame.strip().upper()
    key = _FRAME_ALIASES.get(key, key)
    if key not in SPECTRAL_FRAMES:
        raise ValueError(f"Unknown spectral frame '{frame}'. Must be one of "
                         + ", ".join(SPECTRAL_FRAMES))
    return key


@dataclass(frozen=True)
class SpectralMedium:
    """
    The medium in which the wavelengths of a spectral axis are expressed.

    Parameters
    ----------
    kind : {'vacuum', 'air'}
        Whether wavelengths are vacuum or air wavelengths.
    refraction_method : str, optional
        Which formula for the refractive index of air the wavelengths
        correspond to; one of the methods accepted by
        `~specutils.utils.wcs_utils.refraction_index`. Only meaningful for
        air. Defaults to ``'Morton2000'`` when ``kind`` is ``'air'``.
    temperature : `~astropy.units.Quantity` ['temperature'], optional
        Air temperature. If not given, the standard temperature of the
        refraction formula is assumed.
    pressure : `~astropy.units.Quantity` ['pressure'], optional
        Air pressure. If not given, standard pressure is assumed.
    humidity : float or `~astropy.units.Quantity`, optional
        Relative humidity, as a fraction between 0 and 1 or a percentage
        `~astropy.units.Quantity`. If not given, dry air is assumed.
    co2 : float, optional
        CO2 concentration in ppm. Only used by the ``'Ciddor1996'`` method.

    Examples
    --------
    >>> from specutils.spectra.spectral_frame import SpectralMedium
    >>> SpectralMedium('vacuum')
    SpectralMedium(kind='vacuum')
    >>> SpectralMedium('air', refraction_method='Ciddor1996', co2=400)
    SpectralMedium(kind='air', refraction_method='Ciddor1996', co2=400)
    """
    kind: str
    refraction_method: str = None
    temperature: u.Quantity = None
    pressure: u.Quantity = None
    humidity: object = None
    co2: float = None

    _KINDS = ('vacuum', 'air')

    def __post_init__(self):
        kind = str(self.kind).strip().lower()
        if kind not in self._KINDS:
            raise ValueError(f"SpectralMedium kind must be one of {self._KINDS}, "
                             f"got '{self.kind}'")
        object.__setattr__(self, 'kind', kind)

        conditions = {k: getattr(self, k) for k in
                      ('refraction_method', 'temperature', 'pressure', 'humidity', 'co2')}
        if kind == 'vacuum':
            given = [k for k, v in conditions.items() if v is not None]
            if given:
                raise ValueError("Air conditions cannot be set for a vacuum medium: "
                                 + ", ".join(given))
        else:
            if self.refraction_method is None:
                object.__setattr__(self, 'refraction_method', 'Morton2000')
            if self.temperature is not None:
                temperature = u.Quantity(self.temperature)
                if not temperature.unit.is_equivalent(u.K, equivalencies=u.temperature()):
                    raise u.UnitsError("temperature must have units of temperature")
                object.__setattr__(self, 'temperature', temperature)
            if self.pressure is not None:
                pressure = u.Quantity(self.pressure)
                if not pressure.unit.is_equivalent(u.Pa):
                    raise u.UnitsError("pressure must have units of pressure")
                object.__setattr__(self, 'pressure', pressure)
            if self.humidity is not None:
                humidity = u.Quantity(self.humidity)
                if not humidity.unit.is_equivalent(u.one):
                    raise u.UnitsError("humidity must be dimensionless or a percentage")
                humidity = humidity.to_value(u.one)
                if not 0 <= humidity <= 1:
                    raise ValueError("humidity must be between 0 and 1 (or 0 and 100 percent)")
                object.__setattr__(self, 'humidity', float(humidity))

    @classmethod
    def from_input(cls, value):
        """
        Coerce ``value`` into a `SpectralMedium`.

        Accepts `None`, an existing `SpectralMedium`, the strings ``'vacuum'``
        or ``'air'``, or a dictionary of constructor arguments.
        """
        if value is None or isinstance(value, cls):
            return value
        if isinstance(value, str):
            return cls(value)
        if isinstance(value, dict):
            return cls(**value)
        raise TypeError("medium must be a SpectralMedium, 'vacuum', 'air', "
                        f"or a dict, got {type(value).__name__}")

    @property
    def is_vacuum(self):
        return self.kind == 'vacuum'

    @property
    def is_air(self):
        return self.kind == 'air'

    @property
    def refraction_kwargs(self):
        """
        Keyword arguments describing these air conditions, suitable for
        `~specutils.utils.wcs_utils.refraction_index`,
        `~specutils.utils.wcs_utils.vac_to_air` and
        `~specutils.utils.wcs_utils.air_to_vac`.
        """
        if self.is_vacuum:
            return {}
        return {'method': self.refraction_method, 'co2': self.co2,
                'temperature': self.temperature, 'pressure': self.pressure,
                'humidity': self.humidity}

    def to_dict(self):
        """Plain-dictionary representation, with `None` entries dropped."""
        return {k: v for k, v in asdict(self).items() if v is not None}

    def __str__(self):
        return self.kind

    def __repr__(self):
        args = [f"kind={self.kind!r}"]
        for key, value in self.to_dict().items():
            if key != 'kind':
                args.append(f"{key}={value!r}")
        return f"SpectralMedium({', '.join(args)})"
