from ..ui.widgets.plugin import Plugin
from ..third_party.qtpy.QtWidgets import *
from ..third_party.qtpy.QtCore import *
from ..third_party.qtpy.QtGui import *
from ..core.comms import Dispatch, DispatchHandle
from ..ui.widgets.dialogs import LayerArithmeticDialog
from ..interfaces.managers import layer_manager

from astropy.units import spectral_density, spectral
import logging


class StatisticsPlugin(Plugin):
    name = "Statistics"

    def __init__(self, parent=None):
        super(StatisticsPlugin, self).__init__(parent)

    def setup_ui(self):
        self.layout_vertical.setContentsMargins(11, 11, 11, 11)

        # Setup form layout
        self.layout_form = QFormLayout()
        self.layout_form.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)
        self.layout_form.setFormAlignment(Qt.AlignJustify | Qt.AlignTop)
        self.layout_form.setContentsMargins(1, 1, 1, 12)
        self.layout_form.setSpacing(6)

        self.layout_vertical.addLayout(self.layout_form)

        # Setup labels
        self.label_current_layer = QLabel(self)

        self.layout_form.setWidget(0, QFormLayout.LabelRole, self.label_current_layer)

        self.line_edit_current_layer = QLineEdit(self)
        sizePolicy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.line_edit_current_layer.sizePolicy().hasHeightForWidth())
        self.line_edit_current_layer.setSizePolicy(sizePolicy)
        self.line_edit_current_layer.setStyleSheet("QLineEdit{background: #DDDDDD;}")
        self.line_edit_current_layer.setReadOnly(True)

        self.layout_form.setWidget(0, QFormLayout.FieldRole,
                                   self.line_edit_current_layer)

        # Setup tabs
        self.tab_widget_stats = QTabWidget(self)

        # Setup basic tab
        self.tab_basic = QWidget()
        self.layout_vertical_tab_basic = QVBoxLayout(self.tab_basic)
        self.layout_vertical_tab_basic.setContentsMargins(11, 11, 11, 11)
        self.layout_vertical_tab_basic.setSpacing(6)

        self.layout_form_tab_basic = QFormLayout()
        self.layout_form_tab_basic.setFieldGrowthPolicy(
            QFormLayout.ExpandingFieldsGrow)
        self.layout_form_tab_basic.setFormAlignment(Qt.AlignRight |
                                                    Qt.AlignTop |
                                                    Qt.AlignTrailing)
        self.layout_form_tab_basic.setContentsMargins(1, 1, 1, 1)
        self.layout_form_tab_basic.setSpacing(6)

        self.layout_vertical_tab_basic.addLayout(self.layout_form_tab_basic)

        # Setup basic tab labels
        self.label_mean = QLabel(self.tab_basic)
        self.line_edit_mean = QLineEdit(self.tab_basic)
        self.line_edit_mean.setStyleSheet(
            "QLineEdit{background: #DDDDDD;}")
        self.line_edit_mean.setReadOnly(True)
        self.layout_form_tab_basic.setWidget(0, QFormLayout.LabelRole,
                                             self.label_mean)
        self.layout_form_tab_basic.setWidget(0, QFormLayout.FieldRole,
                                             self.line_edit_mean)

        self.label_median = QLabel(self.tab_basic)
        self.line_edit_median = QLineEdit(self.tab_basic)
        self.line_edit_median.setStyleSheet(
            "QLineEdit{background: #DDDDDD;}")
        self.line_edit_std_dev.setReadOnly(True)
        self.layout_form_tab_basic.setWidget(1, QFormLayout.LabelRole,
                                             self.label_median)
        self.layout_form_tab_basic.setWidget(1, QFormLayout.FieldRole,
                                             self.line_edit_std_dev)

        self.label_std_dev = QLabel(self.tab_basic)
        self.line_edit_std_dev = QLineEdit(self.tab_basic)
        self.line_edit_std_dev.setStyleSheet(
            "QLineEdit{background: #DDDDDD;}")
        self.line_edit_std_dev.setReadOnly(True)
        self.layout_form_tab_basic.setWidget(2, QFormLayout.LabelRole,
                                             self.label_std_dev)
        self.layout_form_tab_basic.setWidget(2, QFormLayout.FieldRole,
                                             self.line_edit_std_dev)

        self.label_total = QLabel(self.tab_basic)
        self.line_edit_total = QLineEdit(self.tab_basic)
        self.line_edit_total.setStyleSheet(
            "QLineEdit{background: #DDDDDD;}")
        self.line_edit_total.setReadOnly(True)
        self.layout_form_tab_basic.setWidget(3, QFormLayout.LabelRole,
                                             self.label_total)
        self.layout_form_tab_basic.setWidget(3, QFormLayout.FieldRole,
                                             self.line_edit_total)

        self.label_data_point_count = QLabel(self.tab_basic)
        self.line_edit_data_point_count = QLineEdit(self.tab_basic)
        self.line_edit_data_point_count.setStyleSheet(
            "QLineEdit{background: #DDDDDD;}")
        self.line_edit_data_point_count.setReadOnly(True)
        self.layout_form_tab_basic.setWidget(4, QFormLayout.LabelRole,
                                             self.label_data_point_count)
        self.layout_form_tab_basic.setWidget(4, QFormLayout.FieldRole,
                                             self.line_edit_data_point_count)


    def setup_connections(self):
        # -- Communications setup
        # Listen for layer selection events, update model tree on selection
        self.tree_widget_layer_list.itemSelectionChanged.connect(
            lambda: Dispatch.on_selected_layer.emit(
                layer_item=self.current_layer_item))

        # When a layer is selected, make that line more obvious than the others
        self.tree_widget_layer_list.itemSelectionChanged.connect(
            lambda: Dispatch.on_selected_plot.emit(
                layer=self.current_layer))

        # When an interactable widget inside a layer item is clicked
        self.tree_widget_layer_list.itemClicked.connect(
            lambda li, col: Dispatch.on_clicked_layer.emit(
                layer_item=li))

        # When an interactable widget inside a layer item is clicked
        self.tree_widget_layer_list.itemChanged.connect(
            lambda li, col: Dispatch.on_changed_layer.emit(
                layer_item=li))

        # -- Widget connection setup
        # When the layer list delete button is pressed
        self.button_remove_layer.clicked.connect(
            lambda: Dispatch.on_remove_layer.emit(layer=self.current_layer))

        # When the arithmetic button is clicked, show math dialog
        self.button_layer_arithmetic.clicked.connect(
            self._show_arithmetic_dialog)

        # Create a new layer based on any active ROIs
        self.button_create_layer_slice.clicked.connect(
            lambda: Dispatch.on_add_roi_layer.emit(layer=self.current_layer, from_roi=True))

    @property
    def current_layer(self):
        """
        Returns the currently selected layer object form the layer list widget.

        Returns
        -------
        layer : specviz.core.data.Layer
            The `Layer` object of the currently selected row.
        """
        layer_item = self.tree_widget_layer_list.currentItem()

        if layer_item is not None:
            layer = layer_item.data(0, Qt.UserRole)

            return layer

    @property
    def current_layer_item(self):
        return self.tree_widget_layer_list.currentItem()

    @property
    def all_layers(self):
        layers = []
        root = self.tree_widget_layer_list.invisibleRootItem()

        for i in range(root.childCount()):
            child = root.child(i)

            if child.data(0, Qt.UserRole):
                layers.append(child.data(0, Qt.UserRole))

            for j in range(child.childCount()):
                sec_child = child.child(j)

                if sec_child.data(0, Qt.UserRole):
                    layers.append(sec_child.data(0, Qt.UserRole))

        return layers

    @DispatchHandle.register_listener("on_added_layer")
    def add_layer_item(self, layer, unique=True):
        """
        Adds a `Layer` object to the loaded layer list widget.

        Parameters
        ----------
        layer : specviz.core.data.Layer
            The `Layer` object to add to the list widget.
        """
        # Make sure there is only one item per layer object
        if unique:
            if self.get_layer_item(layer) is not None:
                return

        new_item = QTreeWidgetItem(
            self.get_layer_item(layer._parent) or
            self.tree_widget_layer_list)
        new_item.setFlags(
            new_item.flags() | Qt.ItemIsUserCheckable | Qt.ItemIsEditable)
        new_item.setText(0, layer.name)
        new_item.setData(0, Qt.UserRole, layer)
        new_item.setCheckState(0, Qt.Checked)

        self.tree_widget_layer_list.setCurrentItem(new_item)

    def get_layer_item(self, layer):
        root = self.tree_widget_layer_list.invisibleRootItem()

        for i in range(root.childCount()):
            child = root.child(i)

            if child.data(0, Qt.UserRole) == layer:
                return child

            for j in range(child.childCount()):
                sec_child = child.child(j)

                if sec_child.data(0, Qt.UserRole) == layer:
                    return sec_child

    @DispatchHandle.register_listener("on_removed_layer")
    def remove_layer_item(self, layer):
        root = self.tree_widget_layer_list.invisibleRootItem()

        for i in range(root.childCount()):
            child = root.child(i)

            if child.data(0, Qt.UserRole) == layer:
                root.removeChild(child)
                break

            for j in range(child.childCount()):
                sec_child = child.child(j)

                if sec_child.data(0, Qt.UserRole) == layer:
                    child.removeChild(sec_child)
                    break

    @DispatchHandle.register_listener("on_added_plot", "on_updated_plot")
    def update_layer_item(self, container=None, *args, **kwargs):
        if container is None:
            return

        layer = container._layer
        pixmap = QPixmap(10, 10)
        pixmap.fill(container.pen.color())
        icon = QIcon(pixmap)

        layer_item = self.get_layer_item(layer)

        if layer_item is not None:
            layer_item.setIcon(0, icon)

    def _show_arithmetic_dialog(self):
        if self.current_layer is None:
            return

        if self.dialog_layer_arithmetic.exec_():
            formula = self.dialog_layer_arithmetic\
                .line_edit_formula.text()

            current_window = self.viewer.current_sub_window
            current_layers = self.all_layers
            new_layer = layer_manager.add_from_formula(formula,
                                                       layers=current_layers)

            if new_layer is None:
                logging.warning("Formula not valid.")
                return

            # If units match, plot the resultant on the same sub window,
            # otherwise create a new sub window to plot the spectra
            data_units_equiv = new_layer.data.unit.is_equivalent(
                current_window._plot_units[1],
                equivalencies=spectral_density(new_layer.dispersion))

            disp_units_equiv = new_layer.dispersion.unit.is_equivalent(
                current_window._plot_units[0], equivalencies=spectral())

            if data_units_equiv and disp_units_equiv:
                self.add_sub_window(layer=new_layer, window=current_window)
            else:
                logging.info("{} not equivalent to {}.".format(
                    new_layer.data.unit, current_window._plot_units[1]))
                self.add_sub_window(layer=new_layer)


