from ..ui.widgets.plugin import Plugin
from ..third_party.qtpy.QtWidgets import *
from ..third_party.qtpy.QtCore import *
from ..third_party.qtpy.QtGui import *
from ..core.comms import Dispatch, DispatchHandle
from ..ui.widgets.dialogs import LayerArithmeticDialog
from ..core.data import Layer

from ..ui.widgets.utils import ICON_PATH

from astropy.units import spectral_density, spectral
import logging


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
        self.layout_vertical_main_widget.setContentsMargins(11, 11, 11, 11)
        self.layout_vertical_main_widget.setSpacing(6)

        self.scroll_area.setWidget(self.main_widget)

        # Tree widget/model selector group box
        self.group_box_add_model = QGroupBox()
        self.group_box_add_model.setTitle("Add Model")
        self.layout_horizontal_group_box_add_model = QHBoxLayout(self.group_box_add_model)
        self.layout_horizontal_group_box_add_model.setContentsMargins(11, 11, 11, 11)
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

        self.layout_horizontal_group_box_add_model.addWidget(self.combo_box_models)
        self.layout_horizontal_group_box_add_model.addWidget(self.button_select_model)

        self.layout_vertical_main_widget.addWidget(self.group_box_add_model)

        # Current models group box
        self.group_box_current_models = QGroupBox(self.main_widget)
        self.group_box_current_models.setTitle("Current Models")
        self.layout_vertical_group_box_current_models = QVBoxLayout(self.group_box_current_models)
        self.layout_vertical_group_box_current_models.setContentsMargins(11, 11,
                                                                         11, 11)
        self.layout_vertical_group_box_current_models.setSpacing(6)

        self.tree_widget_current_models = QTreeWidget(self.group_box_current_models)
        self.tree_widget_current_models.setMinimumSize(QSize(0, 150))
        self.tree_widget_current_models.setAllColumnsShowFocus(False)
        self.tree_widget_current_models.setHeaderHidden(False)
        self.tree_widget_current_models.setColumnCount(2)
        self.tree_widget_current_models.headerItem().setText(0, "Parameter")
        self.tree_widget_current_models.headerItem().setText(1, "Value")
        self.tree_widget_current_models.header().setVisible(True)
        self.tree_widget_current_models.header().setCascadingSectionResizes(False)
        self.tree_widget_current_models.header().setDefaultSectionSize(130)

        self.layout_vertical_group_box_current_models.addWidget(self.tree_widget_current_models)

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
        self.layout_horizontal_model_buttons.addWidget(self.button_export_model)
        self.layout_horizontal_model_buttons.addStretch()
        self.layout_horizontal_model_buttons.addWidget(self.button_remove_model)

        self.layout_vertical_group_box_current_models.addLayout(
            self.layout_horizontal_model_buttons)

        # Arithmetic group box
        self.group_box_model_arithmetic = QGroupBox(self.group_box_current_models)
        self.group_box_model_arithmetic.setTitle("Arithmetic")
        self.layout_vertical_model_arithmetic = QVBoxLayout(self.group_box_model_arithmetic)
        self.layout_vertical_model_arithmetic.setContentsMargins(11, 11, 11,
                                                                 11)
        self.layout_vertical_model_arithmetic.setSpacing(6)

        self.line_edit_model_arithmetic = QLineEdit(self.group_box_model_arithmetic)
        self.layout_vertical_model_arithmetic.addWidget(self.line_edit_model_arithmetic)

        self.layout_vertical_group_box_current_models.addWidget(self.group_box_model_arithmetic)

        # Fitting routines group box
        self.group_box_fitting = QGroupBox(self.main_widget)
        self.group_box_fitting.setTitle("Fitting")
        self.group_box_fitting.setEnabled(False)

        self.layout_vertical_fitting = QVBoxLayout(self.group_box_fitting)
        self.layout_vertical_fitting.setContentsMargins(11, 11, 11, 11)
        self.layout_vertical_fitting.setSpacing(6)

        self.combo_box_fitting = QComboBox(self.group_box_fitting)

        self.button_perform_fit = QPushButton(self.group_box_fitting)

        self.layout_vertical_fitting.addWidget(self.combo_box_fitting)
        self.layout_vertical_fitting.addWidget(self.button_perform_fit)

        # Add group boxees
        self.layout_vertical_main_widget.addWidget(self.group_box_add_model)
        self.layout_vertical_main_widget.addWidget(self.group_box_current_models)
        self.layout_vertical_main_widget.addWidget(
            self.group_box_fitting)

        self.layout_vertical.addWidget(self.scroll_area)

    def setup_connections(self):
        pass