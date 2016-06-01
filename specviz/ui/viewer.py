from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import importlib, inspect, sys, os

from ..third_party.qtpy.QtCore import *
from ..third_party.qtpy.QtWidgets import *
from ..third_party.qtpy.QtGui import *

from .widgets.sub_windows import PlotSubWindow
from ..core.comms import Dispatch, DispatchHandle

from .widgets.windows import MainWindow
from .widgets.plugin import Plugin


class Viewer(object):
    """
    The `Viewer` is the main construction area for all GUI widgets. This
    object does **not** control the interactions between the widgets,
    but only their creation and placement.
    """
    def __init__(self):
        # Instantiate main window object
        self.main_window = MainWindow()

        # Load system and user plugins
        self.load_plugins()

        # Setup up top-level connections
        self._setup_connections()

    def load_plugins(self):
        for mod in os.listdir(os.path.abspath(os.path.join(
                __file__, '..', '..', 'plugins'))):

            mod = mod.split('.')[0]
            mod = importlib.import_module("specviz.plugins." + mod)
            cls_members = inspect.getmembers(
                mod, lambda member: inspect.isclass(member)
                                    and Plugin in member.__bases__)

            if len(cls_members) == 0:
                continue

            for _, cls_plugin in cls_members:
                instance_plugin = cls_plugin(self.main_window)

                if instance_plugin.location != 'hidden':
                    if instance_plugin.location == 'right':
                        location = Qt.RightDockWidgetArea
                    else:
                        location = Qt.LeftDockWidgetArea

                    self.main_window.addDockWidget(location,
                                                   instance_plugin)

    def _setup_connections(self):
        # Listen for subwindow selection events, update layer list on selection
        self.main_window.mdi_area.subWindowActivated.connect(
            lambda wi: Dispatch.on_selected_window.emit(
            window=wi.widget() if wi is not None else None))

        # When a user edits the model parameter field, validate the input
        # self.wgt_model_list.itemChanged.connect(
        #         self._model_parameter_validation)

    @DispatchHandle.register_listener("on_selected_layer")
    def _set_model_tool_options(self, layer_item):
        if layer_item is None:
            self.main_window.createModelLayerButton.hide()
            self.main_window.updateModelLayerButton.hide()
            self.main_window.fittingRoutinesGroupBox.setEnabled(False)
            self.main_window.loadModelButton.setEnabled(False)
            self.main_window.saveModelButton.setEnabled(False)
            self.main_window.exportModelButton.setEnabled(False)

            return

        layer = layer_item.data(0, Qt.UserRole)

        if not hasattr(layer, 'model'):
            self.main_window.createModelLayerButton.show()
            self.main_window.updateModelLayerButton.hide()
            self.main_window.fittingRoutinesGroupBox.setEnabled(False)
            self.main_window.saveModelButton.setEnabled(False)
            self.main_window.exportModelButton.setEnabled(False)
            self.main_window.loadModelButton.setEnabled(True)
        else:
            self.main_window.createModelLayerButton.hide()
            self.main_window.updateModelLayerButton.show()
            self.main_window.fittingRoutinesGroupBox.setEnabled(True)
            self.main_window.saveModelButton.setEnabled(True)
            self.main_window.exportModelButton.setEnabled(True)
            self.main_window.loadModelButton.setEnabled(False)

    @property
    def current_model_item(self):
        return self.wgt_model_list.currentItem()

    @property
    def current_model(self):
        return self.main_window.modelsComboBox.currentText()

    @property
    def current_fitter(self):
        return self.main_window.fittingRoutinesComboBox.currentText()

    @property
    def current_model_formula(self):
        return self.main_window.lineEdit.text()

    @DispatchHandle.register_listener("on_added_model")
    def add_model_item(self, model, layer, unique=True):
        """
        Adds an `astropy.modeling.Model` to the loaded model tree widget.

        Parameters
        ----------
        """
        if model is None:
            return

        if unique:
            if self.get_model_item(model) is not None:
                return

        name = model.name

        if not name:
            count = 1

            root = self.wgt_model_list.invisibleRootItem()

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
            new_para_item.setFlags(new_para_item.flags() | Qt.ItemIsEditable)

        self.wgt_model_list.addTopLevelItem(new_item)

    @DispatchHandle.register_listener("on_removed_model")
    def remove_model_item(self, model=None, layer=None):
        root = self.wgt_model_list.invisibleRootItem()

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

    def get_model_item(self, model):
        root = self.wgt_model_list.invisibleRootItem()

        for i in range(root.childCount()):
            child = root.child(i)

            if child.data(0, Qt.UserRole) == model:
                return child

    def _model_parameter_validation(self, item, col):
        if col == 0:
            return

        try:
            txt = "{:4.4g}".format(float(item.text(col)))
            item.setText(col, txt)
            item.setData(col, Qt.UserRole, float(item.text(col)))
        except ValueError:
            prev_val = item.data(col, Qt.UserRole)
            item.setText(col, str(prev_val))

    def get_model_inputs(self):
        """
        Returns the model and current parameters displayed in the UI.

        Returns
        -------
        models : dict
            A dictionary with the model instance as the key and a list of
            floats as the parameters values.
        """
        root = self.wgt_model_list.invisibleRootItem()
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
