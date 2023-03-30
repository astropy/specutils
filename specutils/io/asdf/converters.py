"""Contains classes that serialize spectral data types into ASDF representations."""
from asdf.extension import Converter
from asdf_astropy.converters import SpectralCoordConverter
from astropy.nddata import (StdDevUncertainty, VarianceUncertainty,
                            InverseVariance, UnknownUncertainty)

from specutils.spectra import Spectrum1D, SpectrumList

__all__ = ['Spectrum1DConverter', 'SpectrumListConverter']

UNCERTAINTY_TYPE_MAPPING = {
    'std': StdDevUncertainty,
    'var': VarianceUncertainty,
    'ivar': InverseVariance,
    'unknown': UnknownUncertainty}


class SpectralAxisConverter(SpectralCoordConverter):
    """ASDF converter to serialize/deserialize SpectralAxis objects."""
    tags = ["tag:astropy.org:specutils/spectra/spectral_axis-*"]
    types = ["specutils.spectra.spectral_axis.SpectralAxis"]

    def from_yaml_tree(self, node, tag, ctx):
        from specutils.spectra.spectral_axis import SpectralAxis

        return SpectralAxis(super().from_yaml_tree(node, tag, ctx))


class Spectrum1DConverter(Converter):
    """ASDF converter to serialize/deserialize Spectrum1D objects."""
    tags = ["tag:astropy.org:specutils/spectra/spectrum1d-*"]
    types = ["specutils.spectra.spectrum1d.Spectrum1D"]

    def to_yaml_tree(self, obj, tag, ctx):
        """Converts Spectrum1D object into tree used for YAML representation."""
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
        """Converts tree representation back into Spectrum1D object."""
        flux = node['flux']
        spectral_axis = node['spectral_axis']
        uncertainty = node.get('uncertainty', None)
        mask = node.get('mask', None)

        if uncertainty is not None:
            class_ = UNCERTAINTY_TYPE_MAPPING[uncertainty['uncertainty_type']]
            data = uncertainty['data']
            uncertainty = class_(data)

        return Spectrum1D(flux=flux, spectral_axis=spectral_axis, uncertainty=uncertainty, mask=mask)


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
