from ..ui.widgets.plugin import Plugin
from ..third_party.qtpy.QtWidgets import *
from ..third_party.qtpy.QtCore import *
from ..third_party.qtpy.QtGui import *
from ..core.comms import Dispatch, DispatchHandle
from ..core.data import GenericSpectrum1DModelLayer
from ..interfaces.model_io import yaml_model_io, py_model_io
from ..interfaces.initializers import initialize
from ..interfaces.factories import ModelFactory, FitterFactory
from ..core.threads import FitModelThread

from ..ui.widgets.utils import ICON_PATH

import numpy as np
import logging

# To memorize last visited directory.
_model_directory = os.environ["HOME"]


class ModelFittingPlugin(Plugin):
    name = "Model Fitting"
    location = "right"

    def __init__(self, *args, **kwargs):
        super(ModelFittingPlugin, self).__init__(*args, **kwargs)
        self.fit_model_thread = FitModelThread()

        self.fit_model_thread.status.connect(
            Dispatch.on_status_message.emit)

        self.fit_model_thread.result.connect(
            lambda layer: Dispatch.on_update_model.emit(layer=layer))

    def setup_ui(self):
        UiModelFittingPlugin(self)

    def setup_connections(self):
        # Enable/disable buttons depending on selection
        self.tree_widget_current_models.itemSelectionChanged.connect(
            self.toggle_buttons)

        # # Populate model dropdown
        self.combo_box_models.addItems(
            sorted(ModelFactory.all_models))

        # Populate fitting algorithm dropdown
        self.combo_box_fitting.addItems(
            sorted(FitterFactory.all_fitters))

        # When the add new model button is clicked, create a new model
        self.button_select_model.clicked.connect(
            self.add_model)

        # When the model items in the model tree change
        self.tree_widget_current_models.itemChanged.connect(
            self._model_parameter_validation)

        # When the model list delete button is pressed
        self.button_remove_model.clicked.connect(
            lambda: self.remove_model_item())

        # When editing the formula is finished, send event
        self.line_edit_model_arithmetic.textEdited.connect(
            lambda: self.update_model_formula())

        # Attach the fit button
        self.button_perform_fit.clicked.connect(
            self.fit_model_layer)

        # ---
        # IO
        # Attach the model save/read buttons
        self.button_save_model.clicked.connect(
            self.save_model)

        self.button_load_model.clicked.connect(
            self.load_model)

        self.button_export_model.clicked.connect(
            self.export_model)

    @property
    def current_model(self):
        model_item = self.current_model_item
        model = model_item.data(0, Qt.UserRole)

        return model

    @property
    def current_model_item(self):
        return self.tree_widget_current_models.currentItem()

    def add_model(self):
        layer = self.current_layer

        if layer is None:
            return

        model_name = self.combo_box_models.currentText()
        model = ModelFactory.create_model(model_name)()

        if isinstance(layer, GenericSpectrum1DModelLayer):
            mask = self.active_window.get_roi_mask(layer._parent)

            initialize(model, layer._parent.dispersion[mask].compressed(),
                       layer._parent.data[mask].compressed())
            # The layer is a `ModelLayer`, in which case, additionally
            # add the model to the compound model and update plot
            layer.model = layer.model + model
        else:
            mask = self.active_window.get_roi_mask(layer)

            initialize(model, layer.dispersion[mask].compressed(),
                       layer.data[mask].compressed())

            # If a layer is selected, but it's not a `ModelLayer`,
            # create a new `ModelLayer`
            layer = self.add_model_layer(model=model)

        Dispatch.on_update_model.emit(layer=layer)
        Dispatch.on_add_model.emit(layer=layer)

    def add_model_layer(self, model):
        """
        Creates a new layer object using the currently defined model.
        """
        layer = self.current_layer

        if layer is None:
            return

        # compound_model = self.get_compound_model(
        #     model_dict=model_inputs,
        #     formula=self.line_edit_model_arithmetic.text())

        # Create new layer using current ROI masks, if they exist
        mask = self.active_window.get_roi_mask(layer=layer)

        new_model_layer = GenericSpectrum1DModelLayer.from_parent(
            parent=layer,
            model=model,
            layer_mask=mask)

        Dispatch.on_add_layer.emit(layer=new_model_layer,
                                   window=self.active_window)

        return new_model_layer

    @DispatchHandle.register_listener("on_add_model")
    def add_model_item(self, layer, unique=True):
        """
        Adds an `astropy.modeling.Model` to the loaded model tree widget.

        Parameters
        ----------
        """
        if hasattr(layer.model, '_submodels'):
            models = layer.model._submodels
        else:
            models = [layer.model]

        for model in models:
            if unique:
                if self.get_model_item(model) is not None:
                    continue

            name = model.name

            if not name:
                count = 1

                root = self.tree_widget_current_models.invisibleRootItem()

                for i in range(root.childCount()):
                    child = root.child(i)

                    if isinstance(model, child.data(0, Qt.UserRole).__class__):
                        count += 1

                name = model.__class__.__name__.replace('1D', '') + str(count)
                model._name = name

            new_item = QTreeWidgetItem()
            new_item.setFlags(new_item.flags() | Qt.ItemIsEditable)

            new_item.setText(0, name)
            new_item.setData(0, Qt.UserRole, model)

            for i, para in enumerate(model.param_names):
                new_para_item = QTreeWidgetItem(new_item)
                new_para_item.setText(0, para)
                new_para_item.setData(1, Qt.UserRole, model.parameters[i])
                new_para_item.setText(1, "{:4.4g}".format(model.parameters[i]))
                new_para_item.setFlags(
                    new_para_item.flags() | Qt.ItemIsEditable)

            self.tree_widget_current_models.addTopLevelItem(new_item)
            self.tree_widget_current_models.expandItem(new_item)

        self._update_arithmetic_text(layer)

    @DispatchHandle.register_listener("on_update_model")
    def update_model_item(self, layer):
        if hasattr(layer.model, '_submodels'):
            models = layer.model._submodels
        else:
            models = [layer.model]

        for model in models:
            model_item = self.get_model_item(model)

            if model_item is None:
                return

            for i, para in enumerate(model.param_names):
                for i in range(model_item.childCount()):
                    param_item = model_item.child(i)

                    if param_item.text(0) == para:
                        param_item.setText(1, "{:4.4g}".format(
                            model.parameters[i]))

    @DispatchHandle.register_listener("on_remove_model")
    def remove_model_item(self, model=None):
        if model is None:
            model = self.current_model

        # Remove model from submodels of compound model
        layer = self.current_layer

        if hasattr(layer, '_model') and hasattr(layer.model, '_submodels'):
            layer.model._submodels.remove(model)
        else:
            logging.error("Cannot remove last model from a `ModelLayer`.")
            return

        # Remove model from tree widget
        root = self.tree_widget_current_models.invisibleRootItem()

        for i in range(root.childCount()):
            child = root.child(i)

            if child.data(0, Qt.UserRole) == model:
                root.removeChild(child)
                break

            for j in range(child.childCount()):
                sec_child = child.child(j)

                if sec_child.data(0, Qt.UserRole) == model:
                    child.removeChild(sec_child)
                    break

        self.update_model_formula()

    def get_model_item(self, model):
        root = self.tree_widget_current_models.invisibleRootItem()

        for i in range(root.childCount()):
            child = root.child(i)

            if child.data(0, Qt.UserRole) == model:
                return child

    def get_model_inputs(self):
        """
        Returns the model and current parameters displayed in the UI.

        Returns
        -------
        models : dict
            A dictionary with the model instance as the key and a list of
            floats as the parameters values.
        """
        root = self.tree_widget_current_models.invisibleRootItem()
        models = {}

        for model_item in [root.child(j) for j in range(root.childCount())]:
            model = model_item.data(0, Qt.UserRole)
            args = []

            for i in range(model_item.childCount()):
                child_item = model_item.child(i)
                child = child_item.text(1)

                args.append(float(child))

            models[model] = args

        return models

    def update_model_formula(self):
        model_layer = self.current_layer

        model_dict = self.get_model_inputs()
        model = self.get_compound_model(model_dict=model_dict)

        if model is None:
            return

        model_layer.model = model

        Dispatch.on_update_model.emit(layer=model_layer)

    def get_compound_model(self, model_dict=None, formula=''):
        model_dict = model_dict or self.get_model_inputs()
        formula = formula or self.line_edit_model_arithmetic.text()
        models = []

        for model in model_dict:
            for i, param_name in enumerate(model.param_names):
                setattr(model, param_name, model_dict[model][i])

            models.append(model)

        if formula:
            model = GenericSpectrum1DModelLayer.from_formula(models, formula)
            return model

        return np.sum(models) if len(models) > 1 else models[0]

    @DispatchHandle.register_listener("on_update_model")
    def _update_arithmetic_text(self, layer):
        if hasattr(layer, '_model'):
            # If the model is a compound
            if hasattr(layer.model, '_submodels'):
                expr = layer.model._format_expression()
                expr = expr.replace('[', '{').replace(']', '}')

                model_names = [model.name
                               for model in layer.model._submodels]

                expr = expr.format(*model_names)
            # If it's just a single model
            else:
                expr = layer.model.name

            self.line_edit_model_arithmetic.setText(expr)

            return expr

    @DispatchHandle.register_listener("on_selected_model", "on_changed_model")
    def _update_model_name(self, model_item, col=0):
        if model_item is None:
            return

        model = model_item.data(0, Qt.UserRole)

        if hasattr(model, '_name'):
            model._name = model_item.text(0)

        self._update_arithmetic_text(self.current_layer)

    @DispatchHandle.register_listener("on_changed_model")
    def _update_model_parameters(self, *args, **kwargs):
        model_layer = self.current_layer
        model_dict = self.get_model_inputs()

        model = self.get_compound_model(model_dict=model_dict,
                                        formula=self.line_edit_model_arithmetic.text())

        if model is not None:
            model_layer.model = model

            Dispatch.on_update_model.emit(layer=model_layer)
        else:
            logging.error("Cannot set `ModelLayer` model to new compound "
                          "model.")

    @DispatchHandle.register_listener("on_selected_layer")
    def update_model_list(self, layer_item=None, layer=None):
        self.tree_widget_current_models.clear()
        self.line_edit_model_arithmetic.clear()

        if layer_item is None and layer is None:
            return

        layer = layer or layer_item.data(0, Qt.UserRole)

        if not hasattr(layer, '_model'):
            return

        self.add_model_item(layer)

    def _model_parameter_validation(self, model_item, col=1):
        if col == 0:
            return

        try:
            txt = "{:4.4g}".format(float(model_item.text(col)))
            model_item.setText(col, txt)
            model_item.setData(col, Qt.UserRole, float(model_item.text(col)))
        except ValueError:
            prev_val = model_item.data(col, Qt.UserRole)
            model_item.setText(col, str(prev_val))

        Dispatch.on_changed_model.emit(model_item=model_item)

    def fit_model_layer(self):
        current_layer = self.current_layer

        if not isinstance(current_layer, GenericSpectrum1DModelLayer):
            logging.error("Attempting to fit model on a non ModelLayer.")
            return

        # This would allow updating the mask on the model layer to reflect
        # the current rois on the plot. Useful for directing fitting,
        # but may be unintuitive.
        mask = self.active_window.get_roi_mask(layer=current_layer._parent)
        current_layer._layer_mask = mask

        # Update the model parameters with those in the gui
        # self.update_model_layer()

        # Create fitted layer
        self.fit_model_thread(
            model_layer=current_layer,
            fitter_name=self.combo_box_fitting.currentText())

        self.fit_model_thread.start()

    def toggle_buttons(self):
        root = self.tree_widget_current_models.invisibleRootItem()

        if root.childCount() > 0:
            self.button_remove_model.setEnabled(True)
        else:
            self.button_remove_model.setEnabled(False)

    @DispatchHandle.register_listener("on_add_model", "on_remove_model")
    def toggle_fitting(self, *args, **kwargs):
        root = self.tree_widget_current_models.invisibleRootItem()

        if root.childCount() > 0:
            self.group_box_fitting.setEnabled(True)
            self.button_save_model.setEnabled(True)
        else:
            self.group_box_fitting.setEnabled(False)
            self.button_save_model.setEnabled(False)

    @DispatchHandle.register_listener("on_selected_layer")
    def toggle_io(self, layer_item, *args, **kwargs):
        if layer_item:
            self.button_load_model.setEnabled(True)
        else:
            self.button_load_model.setEnabled(False)

    # ---
    # IO
    def _prepare_model_for_save(self):
        model_dict = self.get_model_inputs()
        formula = self.line_edit_model_arithmetic.text()

        if len(model_dict) == 0:
            return None, None

        return self.get_compound_model(model_dict,
                                       formula=formula), formula

    def save_model(self):
        model, formula = self._prepare_model_for_save()

        if model:
            global _model_directory
            yaml_model_io.saveModelToFile(self,
                                          model,
                                          _model_directory,
                                          expression=formula)

    def export_model(self):
        model, formula = self._prepare_model_for_save()

        if model:
            global _model_directory
            py_model_io.saveModelToFile(self,
                                        model,
                                        _model_directory,
                                        expression=formula)

    def load_model(self):
        global _model_directory
        fname = QFileDialog.getOpenFileNames(
            self,
            'Read model file',
            _model_directory,
            yaml_model_io.MODEL_FILE_FILTER)

        # File dialog returns a tuple with a list of file names.
        # We get the first name from the first tuple element.
        if len(fname[0]) < 1:
            return
        fname = fname[0][0]

        compound_model, formula, _model_directory = yaml_model_io.buildModelFromFile(
            fname)

        # Put new model in its own sub-layer under current layer.
        current_layer = self.current_layer

        if current_layer is None:
            return

        # Create new model layer using current ROI masks, if they exist
        mask = self.active_window.get_roi_mask(layer=current_layer)

        current_window = self.active_window

        # If there already is a model layer, just edit its model
        if hasattr(current_layer, '_model'):
            current_layer.model = compound_model

            self.update_model_list(layer=current_layer)
            Dispatch.on_update_model.emit(layer=current_layer)
        else:
            new_model_layer = GenericSpectrum1DModelLayer(
                source=current_layer._source,
                mask=mask,
                parent=current_layer,
                name="New Model Layer",
                model=compound_model)

            Dispatch.on_add_layer.emit(layer=new_model_layer,
                                       window=current_window)

        # put formula in text edit widget
        self.line_edit_model_arithmetic.setText(formula)


class UiModelFittingPlugin:
    def __init__(self, plugin):
        # Tree widget/model selector group box
        plugin.group_box_add_model = QGroupBox()
        plugin.group_box_add_model.setTitle("Add Model")
        plugin.layout_horizontal_group_box_add_model = QHBoxLayout(
            plugin.group_box_add_model)
        plugin.layout_horizontal_group_box_add_model.setContentsMargins(11, 11,
                                                                        11, 11)
        plugin.layout_horizontal_group_box_add_model.setSpacing(6)

        # Models combo box
        plugin.combo_box_models = QComboBox(plugin.group_box_add_model)
        plugin.combo_box_models.setStyleSheet("""QComboBox {width: 100px;}
        """)
        size_policy = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred)
        size_policy.setHorizontalStretch(1)
        size_policy.setVerticalStretch(0)
        size_policy.setHeightForWidth(
            plugin.combo_box_models.sizePolicy().hasHeightForWidth())
        plugin.combo_box_models.setSizePolicy(size_policy)

        plugin.button_select_model = QPushButton(plugin.group_box_add_model)
        plugin.button_select_model.setText("Select")

        plugin.layout_horizontal_group_box_add_model.addWidget(
            plugin.combo_box_models)
        plugin.layout_horizontal_group_box_add_model.addWidget(
            plugin.button_select_model)

        plugin.layout_vertical.addWidget(plugin.group_box_add_model)

        # Current models group box
        plugin.group_box_current_models = QGroupBox(plugin.contents)
        plugin.group_box_current_models.setTitle("Current Models")
        plugin.layout_vertical_group_box_current_models = QVBoxLayout(
            plugin.group_box_current_models)
        plugin.layout_vertical_group_box_current_models.setContentsMargins(
            11, 11, 11, 11)
        plugin.layout_vertical_group_box_current_models.setSpacing(6)

        plugin.tree_widget_current_models = QTreeWidget(
            plugin.group_box_current_models)
        # self.tree_widget_current_models.setMinimumSize(QSize(0, 150))
        plugin.tree_widget_current_models.setAllColumnsShowFocus(False)
        plugin.tree_widget_current_models.setHeaderHidden(False)
        plugin.tree_widget_current_models.setColumnCount(2)
        plugin.tree_widget_current_models.headerItem().setText(0, "Parameter")
        plugin.tree_widget_current_models.headerItem().setText(1, "Value")
        plugin.tree_widget_current_models.header().setVisible(True)
        plugin.tree_widget_current_models.header().setCascadingSectionResizes(
            False)
        plugin.tree_widget_current_models.header().setDefaultSectionSize(130)

        plugin.layout_vertical_group_box_current_models.addWidget(
            plugin.tree_widget_current_models)

        # Current models buttons
        plugin.layout_horizontal_model_buttons = QHBoxLayout()
        plugin.layout_horizontal_model_buttons.setContentsMargins(1, 1, 1, 12)
        plugin.layout_horizontal_model_buttons.setSpacing(6)

        plugin.button_save_model = QToolButton(plugin.group_box_current_models)
        plugin.button_save_model.setEnabled(False)
        plugin.button_save_model.setIcon(QIcon(os.path.join(
            ICON_PATH, "Save-48.png")))
        plugin.button_save_model.setIconSize(QSize(25, 25))

        plugin.button_load_model = QToolButton(plugin.group_box_current_models)
        plugin.button_load_model.setEnabled(False)
        plugin.button_load_model.setIcon(QIcon(os.path.join(
            ICON_PATH, "Open Folder-48.png")))
        plugin.button_load_model.setIconSize(QSize(25, 25))

        plugin.button_export_model = QToolButton(plugin.group_box_current_models)
        plugin.button_export_model.setEnabled(False)
        plugin.button_export_model.setIcon(QIcon(os.path.join(
            ICON_PATH, "Export-48.png")))
        plugin.button_export_model.setIconSize(QSize(25, 25))

        plugin.button_remove_model = QToolButton(plugin.group_box_current_models)
        plugin.button_remove_model.setEnabled(False)
        plugin.button_remove_model.setIcon(QIcon(os.path.join(
            ICON_PATH, "Delete-48.png")))
        plugin.button_remove_model.setIconSize(QSize(25, 25))

        plugin.layout_horizontal_model_buttons.addWidget(plugin.button_save_model)
        plugin.layout_horizontal_model_buttons.addWidget(plugin.button_load_model)
        plugin.layout_horizontal_model_buttons.addWidget(
            plugin.button_export_model)
        plugin.layout_horizontal_model_buttons.addStretch()
        plugin.layout_horizontal_model_buttons.addWidget(
            plugin.button_remove_model)

        plugin.layout_vertical_group_box_current_models.addLayout(
            plugin.layout_horizontal_model_buttons)

        # Arithmetic group box
        plugin.group_box_model_arithmetic = QGroupBox(
            plugin.group_box_current_models)
        plugin.group_box_model_arithmetic.setTitle("Arithmetic")
        plugin.layout_vertical_model_arithmetic = QVBoxLayout(
            plugin.group_box_model_arithmetic)
        plugin.layout_vertical_model_arithmetic.setContentsMargins(11, 11, 11,
                                                                   11)
        plugin.layout_vertical_model_arithmetic.setSpacing(6)

        plugin.line_edit_model_arithmetic = QLineEdit(
            plugin.group_box_model_arithmetic)
        plugin.layout_vertical_model_arithmetic.addWidget(
            plugin.line_edit_model_arithmetic)

        plugin.layout_vertical_group_box_current_models.addWidget(
            plugin.group_box_model_arithmetic)

        # Fitting routines group box
        plugin.group_box_fitting = QGroupBox(plugin.contents)
        plugin.group_box_fitting.setTitle("Fitting")
        plugin.group_box_fitting.setEnabled(False)

        plugin.layout_vertical_fitting = QVBoxLayout(plugin.group_box_fitting)
        plugin.layout_vertical_fitting.setContentsMargins(11, 11, 11, 11)
        plugin.layout_vertical_fitting.setSpacing(6)

        plugin.combo_box_fitting = QComboBox(plugin.group_box_fitting)

        plugin.button_perform_fit = QPushButton(plugin.group_box_fitting)
        plugin.button_perform_fit.setText("Perform fit")

        plugin.layout_vertical_fitting.addWidget(plugin.combo_box_fitting)
        plugin.layout_vertical_fitting.addWidget(plugin.button_perform_fit)

        # Add group boxees
        plugin.layout_vertical.addWidget(plugin.group_box_add_model)
        plugin.layout_vertical.addWidget(
            plugin.group_box_current_models)
        plugin.layout_vertical.addWidget(
            plugin.group_box_fitting)




