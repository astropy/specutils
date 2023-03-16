"""
Contains classes that serialize spectral data types into ASDF representations.
"""
from asdf.extension import Converter
from asdf.yamlutil import custom_tree_to_tagged_tree, tagged_tree_to_custom_tree
from astropy.nddata import (StdDevUncertainty, VarianceUncertainty,
                            InverseVariance, UnknownUncertainty)

from specutils.spectra import Spectrum1D, SpectrumList

__all__ = ['Spectrum1DType', 'SpectrumListType']

UNCERTAINTY_TYPE_MAPPING = {
    'std': StdDevUncertainty,
    'var': VarianceUncertainty,
    'ivar': InverseVariance,
    'unknown': UnknownUncertainty}


class Spectrum1DType(Converter):
    """ASDF tag implementation used to serialize/deserialize Spectrum1D objects."""
    tags = ["tag:astropy.org:astropy/spectra/spectrum1d-*"]
    types = ["specutils.spectra.spectrum1d.Spectrum1D"]

    def to_yaml_tree(self, obj, tag, ctx):
        """Converts Spectrum1D object into tree used for YAML representation."""
        node = {}
        node['flux'] = custom_tree_to_tagged_tree(obj.flux, ctx)
        node['spectral_axis'] = custom_tree_to_tagged_tree(obj.spectral_axis, ctx)

        if obj.uncertainty is not None:
            node['uncertainty'] = {}
            node['uncertainty']['uncertainty_type'] = obj.uncertainty.uncertainty_type
            data = custom_tree_to_tagged_tree(obj.uncertainty.array, ctx)
            node['uncertainty']['data'] = data

        return node

    def from_yaml_tree(cls, node, tag, ctx):
        """Converts tree representation back into Spectrum1D object."""
        flux = tagged_tree_to_custom_tree(node['flux'], ctx)
        spectral_axis = tagged_tree_to_custom_tree(node['spectral_axis'], ctx)
        uncertainty = node.get('uncertainty', None)
        if uncertainty is not None:
            klass = UNCERTAINTY_TYPE_MAPPING[uncertainty['uncertainty_type']]
            data = tagged_tree_to_custom_tree(uncertainty['data'], ctx)
            uncertainty = klass(data)

        return Spectrum1D(flux=flux, spectral_axis=spectral_axis, uncertainty=uncertainty)


class SpectrumListType(Converter):
    """ASDF tag implementation used to serialize/deserialize SpectrumList objects."""
    tags = ["tag:astropy.org:astropy/spectra/spectrum_list-*"]
    types = ["specutils.spectra.spectrum_list.SpectrumList"]

    def to_yaml_tree(self, obj, tag, ctx):
        """Converts SpectrumList object into tree used for YAML representation."""
        return [custom_tree_to_tagged_tree(spectrum, ctx) for spectrum in obj]

    def from_yaml_tree(cls, node, tag, ctx):
        """Converts tree representation back into SpectrumList object."""
        return SpectrumList(tagged_tree_to_custom_tree(tree, ctx) for tree in node)
