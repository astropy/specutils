import logging

import numpy as np
from astropy.units import spectral_density

from ..analysis.utils import resample


class LayerArithmeticMixin:
    def _arithmetic(self, operator, operand, propagate_uncertainties=True,
                    **kwargs):
        # Make sure the classes match
        operand = self.__class__(operand)

        # Make sure units are compatible
        if not operand.unit.is_equivalent(self.unit,
                                          equivalencies=spectral_density(
                                            operand.dispersion)):
            logging.error("Spectral data objects have incompatible units.")
            return

        # Always be sure that the wavelength grid encompasses both layers
        if operand.dispersion.data.to(self.unit)[-1] > \
                self.dispersion.data[-1]:
            max_disp = operand.dispersion.data.to(self.unit)[-1]
            logging.info("Dispersion grids require extrapolation.")
        else:
            max_disp = self.dispersion.data[-1]

        if operand.dispersion.data.to(self.unit)[0] < \
                self.dispersion.data[0]:
            min_disp = operand.dispersion.data.to(self.unit)[0]
            logging.info("Dispersion grids require extrapolation.")
        else:
            min_disp = self.dispersion.data[0]

        step = self.dispersion.data[1] - self.dispersion.data[0]
        new_dispersion = np.arange(min_disp, max_disp, step)

        # Resample onto new dispersion grid
        other_samp = np.ones(new_dispersion.shape,
                             [('wlen', float), ('flux', float),
                              ('ivar', float)])
        other_samp['wlen'] = operand.dispersion.data.value
        other_samp['flux'] = operand.data.data.value
        other_samp['ivar'] = operand.uncertainty if operand.uncertainty \
                             is not None else np.zeros(
            operand.dispersion.data.value.shape)

        self_samp = np.ones(new_dispersion.shape,
                             [('wlen', float), ('flux', float),
                              ('ivar', float)])
        self_samp['wlen'] = self.dispersion.data.value
        self_samp['flux'] = self.data.data.value
        self_samp['ivar'] = self.uncertainty if operand.uncertainty \
                            is not None else np.zeros(
            self.dispersion.data.value.shape)

        # Define the columns to resample
        cols = ('flux', 'ivar')

        other_resamp = resample(other_samp, 'wlen', new_dispersion, cols)
        self_resamp = resample(self_samp, 'wlen', new_dispersion, cols)

        new_self = self._source.copy(self_resamp, dispersion=new_dispersion)
        new_operand = operand._source.copy(other_resamp,
                                           dispersion=new_dispersion)

        # Perform arithmetic operation
        operator = getattr(new_self, operator)
        result_source = operator(new_operand, propagate_uncertainties=propagate_uncertainties)

        # Generate the final mask
        final_mask = self.mask | operand.mask

        # Create a layer from the source data object
        result_layer = Layer(result_source, mask=~final_mask,
                             parent=self._parent, name=self.name)

        return result_layer

    def __add__(self, other):
        new_layer = self._arithmetic("add", other)

        return new_layer

    def __sub__(self, other):
        new_layer = self._arithmetic("subtract", other)

        return new_layer

    def __mul__(self, other):
        new_layer = self._arithmetic("multiply", other, propagate_uncertainties=True)

        return new_layer

    def __truediv__(self, other):
        new_layer = self._arithmetic("divide", other, propagate_uncertainties=True)

        return new_layer


class ModelLayerArithmeticMixin:
    def _arithmetic(self):
        # In the case where the model layers are from two separate data
        # layers, create a combined data layer
        # operator_char = dict(add='+', subtract='-', divide='-',
        #                      multiply='*')[operator]
        # final_layer = Layer.from_formula(
        #     "{} {} {}".format(self._parent.name, operator_char,
        #                       other._parent.name),
        #     [self._parent, other._parent])
        operator_func = dict(add='add', subtract='sub',
                             divide='truediv', multiply='mult')[operator]
        final_model = getattr(self.model,
                              "__{}__".format(operator_func))(other.model)

        final_names = []

        for smod in final_model._submodels:
            if smod._name in final_names:
                split_name = re.findall(r"[^\W\d_]+|\d+", smod._name)
                split_name[-1] = str(int(split_name[-1]) + 1)
                smod._name = "".join(split_name)

            final_names.append(smod._name)

        final_mask = self.mask | other.full_mask

        result_model_layer = ModelLayer(model=final_model,
                                        source=self._parent._source,
                                        mask=~final_mask,
                                        parent=self._parent)

        return result_model_layer
