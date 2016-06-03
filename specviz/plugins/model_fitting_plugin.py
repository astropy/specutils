from ..ui.widgets.plugin import Plugin
from ..third_party.qtpy.QtWidgets import *
from ..third_party.qtpy.QtCore import *
from ..third_party.qtpy.QtGui import *
from ..core.comms import Dispatch, DispatchHandle
from ..core.data import ModelLayer
# from ..interfaces.managers import model_manager
from ..interfaces.model_io import yaml_model_io, py_model_io
from ..interfaces.factories import ModelFactory, FitterFactory

from ..ui.widgets.utils import ICON_PATH

from astropy.units import spectral_density, spectral
import logging
import numpy as np


class ModelFittingPlugin(Plugin):
    name = "Model Fitting"
    location = "right"

    def setup_ui(self):
        self.scroll_area = QScrollArea(self)
        self.scroll_area.setFrameShape(QFrame.NoFrame)
        self.scroll_area.setFrameShadow(QFrame.Plain)
        self.scroll_area.setLineWidth(0)
        self.scroll_area.setWidgetResizable(True)

        # The main widget inside the scroll area
        self.main_widget = QWidget()
        self.layout_vertical_main_widget = QVBoxLayout(self.main_widget)
        self.layout_vertical_main_widget.setContentsMargins(1, 1, 1, 1)
        self.layout_vertical_main_widget.setSpacing(6)

        self.scroll_area.setWidget(self.main_widget)

        # Tree widget/model selector group box
        self.group_box_add_model = QGroupBox()
        self.group_box_add_model.setTitle("Add Model")
        self.layout_horizontal_group_box_add_model = QHBoxLayout(
            self.group_box_add_model)
        self.layout_horizontal_group_box_add_model.setContentsMargins(11, 11,
                                                                      11, 11)
        self.layout_horizontal_group_box_add_model.setSpacing(6)

        # Models combo box
        self.combo_box_models = QComboBox(self.group_box_add_model)

        size_policy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        size_policy.setHorizontalStretch(1)
        size_policy.setVerticalStretch(0)
        size_policy.setHeightForWidth(
            self.combo_box_models.sizePolicy().hasHeightForWidth())
        self.combo_box_models.setSizePolicy(size_policy)

        self.button_select_model = QPushButton(self.group_box_add_model)
        self.button_select_model.setText("Select")

        self.layout_horizontal_group_box_add_model.addWidget(
            self.combo_box_models)
        self.layout_horizontal_group_box_add_model.addWidget(
            self.button_select_model)

        self.layout_vertical_main_widget.addWidget(self.group_box_add_model)

        # Current models group box
        self.group_box_current_models = QGroupBox(self.main_widget)
        self.group_box_current_models.setTitle("Current Models")
        self.layout_vertical_group_box_current_models = QVBoxLayout(
            self.group_box_current_models)
        self.layout_vertical_group_box_current_models.setContentsMargins(11,
                                                                         11,
                                                                         11,
                                                                         11)
        self.layout_vertical_group_box_current_models.setSpacing(6)

        self.tree_widget_current_models = QTreeWidget(
            self.group_box_current_models)
        self.tree_widget_current_models.setMinimumSize(QSize(0, 150))
        self.tree_widget_current_models.setAllColumnsShowFocus(False)
        self.tree_widget_current_models.setHeaderHidden(False)
        self.tree_widget_current_models.setColumnCount(2)
        self.tree_widget_current_models.headerItem().setText(0, "Parameter")
        self.tree_widget_current_models.headerItem().setText(1, "Value")
        self.tree_widget_current_models.header().setVisible(True)
        self.tree_widget_current_models.header().setCascadingSectionResizes(
            False)
        self.tree_widget_current_models.header().setDefaultSectionSize(130)

        self.layout_vertical_group_box_current_models.addWidget(
            self.tree_widget_current_models)

        # Current models buttons
        self.layout_horizontal_model_buttons = QHBoxLayout()
        self.layout_horizontal_model_buttons.setContentsMargins(1, 1, 1, 12)
        self.layout_horizontal_model_buttons.setSpacing(6)

        self.button_save_model = QToolButton(self.group_box_current_models)
        self.button_save_model.setEnabled(False)
        self.button_save_model.setIcon(QIcon(os.path.join(
            ICON_PATH, "Save-48.png")))
        self.button_save_model.setIconSize(QSize(25, 25))

        self.button_load_model = QToolButton(self.group_box_current_models)
        self.button_load_model.setEnabled(False)
        self.button_load_model.setIcon(QIcon(os.path.join(
            ICON_PATH, "Open Folder-48.png")))
        self.button_load_model.setIconSize(QSize(25, 25))

        self.button_export_model = QToolButton(self.group_box_current_models)
        self.button_export_model.setEnabled(False)
        self.button_export_model.setIcon(QIcon(os.path.join(
            ICON_PATH, "Export-48.png")))
        self.button_export_model.setIconSize(QSize(25, 25))

        self.button_remove_model = QToolButton(self.group_box_current_models)
        self.button_remove_model.setEnabled(False)
        self.button_remove_model.setIcon(QIcon(os.path.join(
            ICON_PATH, "Delete-48.png")))
        self.button_remove_model.setIconSize(QSize(25, 25))

        self.layout_horizontal_model_buttons.addWidget(self.button_save_model)
        self.layout_horizontal_model_buttons.addWidget(self.button_load_model)
        self.layout_horizontal_model_buttons.addWidget(
            self.button_export_model)
        self.layout_horizontal_model_buttons.addStretch()
        self.layout_horizontal_model_buttons.addWidget(
            self.button_remove_model)

        self.layout_vertical_group_box_current_models.addLayout(
            self.layout_horizontal_model_buttons)

        # Arithmetic group box
        self.group_box_model_arithmetic = QGroupBox(
            self.group_box_current_models)
        self.group_box_model_arithmetic.setTitle("Arithmetic")
        self.layout_vertical_model_arithmetic = QVBoxLayout(
            self.group_box_model_arithmetic)
        self.layout_vertical_model_arithmetic.setContentsMargins(11, 11, 11,
                                                                 11)
        self.layout_vertical_model_arithmetic.setSpacing(6)

        self.line_edit_model_arithmetic = QLineEdit(
            self.group_box_model_arithmetic)
        self.layout_vertical_model_arithmetic.addWidget(
            self.line_edit_model_arithmetic)

        self.layout_vertical_group_box_current_models.addWidget(
            self.group_box_model_arithmetic)

        # Fitting routines group box
        self.group_box_fitting = QGroupBox(self.main_widget)
        self.group_box_fitting.setTitle("Fitting")
        self.group_box_fitting.setEnabled(False)

        self.layout_vertical_fitting = QVBoxLayout(self.group_box_fitting)
        self.layout_vertical_fitting.setContentsMargins(11, 11, 11, 11)
        self.layout_vertical_fitting.setSpacing(6)

        self.combo_box_fitting = QComboBox(self.group_box_fitting)

        self.button_perform_fit = QPushButton(self.group_box_fitting)
        self.button_perform_fit.setText("Perform fit")

        self.layout_vertical_fitting.addWidget(self.combo_box_fitting)
        self.layout_vertical_fitting.addWidget(self.button_perform_fit)

        # Add group boxees
        self.layout_vertical_main_widget.addWidget(self.group_box_add_model)
        self.layout_vertical_main_widget.addWidget(
            self.group_box_current_models)
        self.layout_vertical_main_widget.addWidget(
            self.group_box_fitting)

        self.layout_vertical.addWidget(self.scroll_area)

    def setup_connections(self):
        # Populate model dropdown
        self.combo_box_models.addItems(
            ModelFactory.all_models)

        # Populate fitting algorithm dropdown
        self.combo_box_fitting.addItems(
            FitterFactory.all_fitters)

        # When the add new model button is clicked, create a new model
        self.button_select_model.clicked.connect(
            self.add_model)

        # When the model items in the model tree change
        self.tree_widget_current_models.itemChanged.connect(
            lambda mi, col: Dispatch.on_changed_model.emit(
                model_item=mi))

        # When the model list delete button is pressed
        self.button_remove_model.clicked.connect(
            lambda: self.remove_model())

        # When editing the formula is finished, send event
        self.line_edit_model_arithmetic.textEdited.connect(
            lambda: self.update_model_formula())

    @property
    def current_model(self):
        model_item = self.current_model_item()
        model = model_item.data(0, Qt.UserRole)

        return model

    @property
    def current_model_item(self):
        return self.tree_widget_current_models.currentItem()

    def add_model(self):
        model_name = self.combo_box_models.currentText()
        model = ModelFactory.create_model(model_name)()
        layer = self.current_layer

        if layer is not None:
            # There is current a selected layer
            if hasattr(layer, '_model'):
                # The layer is a `ModelLayer`, in which case, additionally
                # add the model to the compound model and update plot
                layer.model = layer.model + model

                mask = self.active_window.get_roi_mask(layer._parent)
            else:
                # If a layer is selected, but it's not a `ModelLayer`,
                # create a new `ModelLayer`
                layer = self.add_model_layer(model=model)
                mask = self.active_window.get_roi_mask(layer)
        else:
            return
            # There is no currently selected layer, simply

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

        new_model_layer = ModelLayer(
            source=layer._source,
            mask=mask,
            parent=layer,
            name="New Model Layer",
            model=model)

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
                new_para_item.setData(0, Qt.UserRole,
                                      model.parameters[i])
                new_para_item.setText(1, "{:4.4g}".format(model.parameters[i]))
                new_para_item.setFlags(
                    new_para_item.flags() | Qt.ItemIsEditable)

            self.tree_widget_current_models.addTopLevelItem(new_item)

        self._update_arithmetic_text(layer)

    @DispatchHandle.register_listener("on_removed_model")
    def remove_model_item(self, model=None, layer=None):
        root = self.tree_widget_current_models.invisibleRootItem()

        for i in range(root.childCount()):
            child = root.child(i)

            if child is None:
                continue

            if child.data(0, Qt.UserRole) == model:
                root.removeChild(child)
                break

    def update_model_item(self, model):
        if hasattr(model, '_submodels'):
            for sub_model in model._submodels:
                self.update_model_item(sub_model)
            else:
                return

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

    def get_compound_model(self, model_dict, formula=''):
        formula = formula or self.line_edit_model_arithmetic.text()
        models = []

        for model in model_dict:
            for i, param_name in enumerate(model.param_names):
                setattr(model, param_name, model_dict[model][i])

            models.append(model)

        if formula:
            model = ModelLayer.from_formula(models, formula)

            return model

        return np.sum(models) if len(models) > 1 else models[0]

    def _update_arithmetic_text(self, model_layer):
        if hasattr(model_layer, '_model'):
            # If the model is a compound
            if hasattr(model_layer.model, '_submodels'):
                expr = model_layer.model._format_expression()
                expr = expr.replace('[', '{').replace(']', '}')

                model_names = [model.name
                               for model in model_layer.model._submodels]

                expr = expr.format(*model_names)
            # If it's just a single model
            else:
                expr = model_layer.model.name

            self.line_edit_model_arithmetic.setText(expr)

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

        model_layer.model = model

        Dispatch.on_update_model.emit(layer=model_layer)

    @DispatchHandle.register_listener("on_selected_layer")
    def update_model_list(self, layer_item):
        self.tree_widget_current_models.clear()
        self.line_edit_model_arithmetic.clear()

        if layer_item is None:
            return

        layer = layer_item.data(0, Qt.UserRole)

        if not hasattr(layer, '_model'):
            return

        self.add_model_item(layer)

    @DispatchHandle.register_listener("on_changed_model")
    def _model_parameter_validation(self, model_item, col):
        if col == 0:
            return

        try:
            txt = "{:4.4g}".format(float(model_item.text(col)))
            model_item.setText(col, txt)
            model_item.setData(col, Qt.UserRole, float(model_item.text(col)))
        except ValueError:
            prev_val = model_item.data(col, Qt.UserRole)
            model_item.setText(col, str(prev_val))

