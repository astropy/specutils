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
