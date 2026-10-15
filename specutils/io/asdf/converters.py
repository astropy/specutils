"""Contains classes that serialize spectral data types into ASDF representations."""
import numpy as np
from asdf.extension import Converter
from asdf_astropy.converters import SpectralCoordConverter
from astropy.coordinates import SkyCoord
from astropy.nddata import (StdDevUncertainty, VarianceUncertainty,
                            InverseVariance, UnknownUncertainty)

from specutils.spectra import Spectrum, SpectrumList

__all__ = ['SpectrumConverter', 'SpectrumListConverter']

UNCERTAINTY_TYPE_MAPPING = {
    'std': StdDevUncertainty,
    'var': VarianceUncertainty,
    'ivar': InverseVariance,
    'unknown': UnknownUncertainty}


class SpectralAxisConverter(SpectralCoordConverter):
    """ASDF converter to serialize/deserialize SpectralAxis objects."""
    tags = ["tag:astropy.org:specutils/spectra/spectral_axis-*"]
    types = ["specutils.spectra.spectral_axis.SpectralAxis"]

    def to_yaml_tree(self, obj, tag, ctx):
        if tag.endswith("-1.0.0"):
            return super().to_yaml_tree(obj, tag, ctx)

        node = {"value": obj.value, "unit": obj.unit}

        derived = obj.observer is not None and obj.target is not None
        if obj._radial_velocity is not None and not derived:
            node["radial_velocity"] = obj.radial_velocity
        if obj.doppler_rest is not None:
            node["doppler_rest"] = obj.doppler_rest
        if obj.doppler_convention is not None:
            node["doppler_convention"] = obj.doppler_convention

        if obj.medium is not None:
            node["medium"] = obj.medium.to_dict()
        if obj.frame is not None:
            node["frame"] = obj.frame
        if obj.obstime is not None:
            node["obstime"] = obj.obstime
        if obj.location is not None:
            node["location"] = obj.location

        # Coordinate frames lose their velocities in ASDF but SkyCoords keep
        # them, so store both as SkyCoord. The observer is only stored if it
        # cannot be rebuilt from the frame, location and time on read.
        if obj.target is not None:
            node["target"] = SkyCoord(obj.target)
        if obj.observer is not None and not (obj.frame is not None
                                             and obj.location is not None
                                             and obj.obstime is not None):
            node["observer"] = SkyCoord(obj.observer)

        return node

    def from_yaml_tree(self, node, tag, ctx):
        from specutils.spectra.spectral_axis import SpectralAxis
        from specutils.spectra.spectral_frame import SpectralMedium

        if tag.endswith("-1.0.0"):
            return SpectralAxis(super().from_yaml_tree(node, tag, ctx))

        medium = node.get("medium")
        if medium is not None:
            medium = SpectralMedium(**medium)

        return SpectralAxis(
            np.asarray(node["value"]), unit=node["unit"],
            radial_velocity=node.get("radial_velocity"),
            doppler_rest=node.get("doppler_rest"),
            doppler_convention=node.get("doppler_convention"),
            medium=medium, frame=node.get("frame"),
            obstime=node.get("obstime"), location=node.get("location"),
            target=node.get("target"), observer=node.get("observer"))


class SpectrumConverter(Converter):
    """ASDF converter to serialize/deserialize Spectrum objects."""
    tags = ["tag:astropy.org:specutils/spectra/spectrum-*",
            "tag:astropy.org:specutils/spectra/spectrum1d-*"]
    types = ["specutils.spectra.spectrum.Spectrum",
             "specutils.spectra.spectrum.Spectrum1D"]

    def to_yaml_tree(self, obj, tag, ctx):
        """Converts Spectrum object into tree used for YAML representation."""
        node = {}
        node['flux'] = obj.flux
        node['spectral_axis'] = obj.spectral_axis

        if obj.uncertainty is not None:
            node['uncertainty'] = {}
            node['uncertainty']['uncertainty_type'] = obj.uncertainty.uncertainty_type
            data = obj.uncertainty.array
            node['uncertainty']['data'] = data

        if obj.mask is not None:
            node['mask'] = obj.mask

        return node

    def from_yaml_tree(cls, node, tag, ctx):
        """Converts tree representation back into Spectrum object."""
        flux = node['flux']
        spectral_axis = node['spectral_axis']
        uncertainty = node.get('uncertainty', None)
        mask = node.get('mask', None)

        if uncertainty is not None:
            class_ = UNCERTAINTY_TYPE_MAPPING[uncertainty['uncertainty_type']]
            data = uncertainty['data']
            uncertainty = class_(data)

        return Spectrum(flux=flux, spectral_axis=spectral_axis, uncertainty=uncertainty, mask=mask)


class SpectrumListConverter(Converter):
    """ASDF converter used to serialize/deserialize SpectrumList objects."""
    tags = ["tag:astropy.org:specutils/spectra/spectrum_list-*"]
    types = ["specutils.spectra.spectrum_list.SpectrumList"]

    def to_yaml_tree(self, obj, tag, ctx):
        """Converts SpectrumList object into tree used for YAML representation."""
        return [spectrum for spectrum in obj]

    def from_yaml_tree(cls, node, tag, ctx):
        """Converts tree representation back into SpectrumList object."""
        return SpectrumList(tree for tree in node)
