from ..third_party.qtpy.QtCore import QThread, pyqtSignal
import os
import logging

from ..core.data import Spectrum1DRef, Spectrum1DRefModelLayer
from ..interfaces.factories import FitterFactory

import astropy.io.registry as io_registry


class FileLoadThread(QThread):
    status = pyqtSignal(str, int)
    result = pyqtSignal(Spectrum1DRef)

    def __init__(self, parent=None):
        super(FileLoadThread, self).__init__(parent)
        self.file_name = ""
        self.file_filter = ""

    def __call__(self, file_name, file_filter):
        self.file_name = file_name
        self.file_filter = file_filter

    def run(self):
        self.status.emit("Loading file...", 0)
        data = self.read_file(self.file_name, self.file_filter)

        if data is not None:
            self.status.emit("File loaded successfully!", 5000)
        else:
            self.status.emit("An error occurred while loading file.", 5000)

        if data is not None:
            self.result.emit(data)
        else:
            logging.error("Could not open file.")

    def read_file(self, file_name, file_filter):
        """
        Convenience method that directly reads a spectrum from a file.
        This exists mostly to facilitate development workflow. In time it
        could be augmented to support fancier features such as wildcards,
        file lists, mixed file types, and the like.
        Note that the filter string is hard coded here; its details might
        depend on the intrincacies of the registries, loaders, and data
        classes. In other words, this is brittle code.
        """
        file_name = str(file_name)
        file_ext = os.path.splitext(file_name)[-1]

        if file_filter == 'Auto':
            all_formats = io_registry.get_formats(Spectrum1DRef)['Format']
        else:
            all_formats = [file_filter]

        for format in all_formats:
            try:
                data = Spectrum1DRef.read(file_name, format=format)
                return data
            except:
                logging.error("Incompatible loader for selected data: {"
                              "}".format(file_filter))


class FitModelThread(QThread):
    status = pyqtSignal(str, int)
    result = pyqtSignal(Spectrum1DRefModelLayer)

    def __init__(self, parent=None):
        super(FitModelThread, self).__init__(parent)
        self.model_layer = None
        self.fitter_name = ""

    def __call__(self, model_layer, fitter_name):
        self.model_layer = model_layer
        self.fitter_name = fitter_name

    def run(self):
        self.status.emit("Fitting model...", 0)
        model_layer, message = self.fit_model(self.model_layer, self.fitter_name)

        if not message:
            self.status.emit("Fit completed successfully!", 5000)
        else:
            self.status.emit("Fit completed, but with warnings.", 5000)

        self.result.emit(model_layer)

    def fit_model(self, model_layer, fitter_name):
        if not hasattr(model_layer, 'model'):
            logging.warning("This layer has no model to fit.")
            return

        # When fitting, the selected layer is a ModelLayer, thus
        # the data to be fitted resides in the parent
        parent_layer = model_layer._parent

        if parent_layer is None:
            return

        flux = parent_layer.data
        dispersion = parent_layer.dispersion
        model = model_layer.model

        # The fitting should only consider the masked regions
        flux = flux[model_layer.layer_mask].compressed().value
        dispersion = dispersion[model_layer.layer_mask].compressed().value

        # Get compressed versions of the data arrays
        # flux = flux.compressed().value
        # dispersion = dispersion.compressed().value

        # If the number of parameters is greater than the number of data
        # points, bail
        if len(model.parameters) > flux.size:
            logging.warning("Unable to perform fit; number of parameters is "
                            "greater than the number of data points.")
            return

        # Perform fitting of model
        if fitter_name:
            fitter = FitterFactory.all_fitters[fitter_name]()
        else:
            fitter = FitterFactory.default_fitter()

        fitted_model = fitter(model, dispersion, flux, maxiter=2000)

        if 'message' in fitter.fit_info:
            # The fitter 'message' should probably be logged at INFO level.
            # Problem is, info messages do not display in the error console,
            # and we, ideally, want the user to see the message immediately
            # after the fit is executed.
            logging.warning(fitter.fit_info['message'])

        # Update original model with new values from fitted model
        if hasattr(fitted_model, '_submodels'):
            for i in range(len(fitted_model._submodels)):
                for pname in model._submodels[i].param_names:
                    value = getattr(fitted_model, "{}_{}".format(pname, i))
                    setattr(model._submodels[i], pname, value.value)
                    setattr(model[i], pname, value.value)
        else:
            for pname in model.param_names:
                value = getattr(fitted_model, "{}".format(pname))
                setattr(model, pname, value.value)
        # model_layer.model = fitted_model

        # update GUI with fit results

        return model_layer, fitter.fit_info.get('message', "")

