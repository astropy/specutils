from ..ui.widgets.plugin import Plugin
from ..third_party.qtpy.QtWidgets import *
from ..third_party.qtpy.QtCore import *
from ..third_party.qtpy.QtGui import *
from ..core.comms import Dispatch, DispatchHandle
from ..ui.widgets.dialogs import LayerArithmeticDialog
from ..interfaces.factories import DataFactory
from ..core.data import Layer

from ..ui.widgets.utils import ICON_PATH

from astropy.units import spectral_density, spectral
import logging


class LayerListPlugin(Plugin):
    name = "Layer List"

    def setup_ui(self):
        self.layout_vertical.setContentsMargins(11, 11, 11, 11)

        self.tree_widget_layer_list = QTreeWidget(self)

        self.layout_vertical.addWidget(self.tree_widget_layer_list)

        self.layout_horizontal = QHBoxLayout()

        self.button_layer_arithmetic = QToolButton(self)
        self.button_layer_arithmetic.setIcon(QIcon(os.path.join(
            ICON_PATH, "Math-48.png")))
        self.button_layer_arithmetic.setEnabled(False)
        self.button_layer_arithmetic.setIconSize(QSize(25, 25))

        self.button_remove_layer = QToolButton(self)
        self.button_remove_layer.setIcon(QIcon(os.path.join(
            ICON_PATH, "Delete-48.png")))
        self.button_remove_layer.setEnabled(False)
        self.button_remove_layer.setIconSize(QSize(25, 25))

        self.button_change_color = QToolButton(self)
        self.button_change_color.setIcon(QIcon(os.path.join(
            ICON_PATH, "Color Dropper-48.png")))
        self.button_change_color.setEnabled(False)
        self.button_change_color.setIconSize(QSize(25, 25))

        self.layout_horizontal.addWidget(self.button_layer_arithmetic)
        self.layout_horizontal.addStretch()
        self.layout_horizontal.addWidget(self.button_change_color)
        self.layout_horizontal.addWidget(self.button_remove_layer)

        self.layout_vertical.addLayout(self.layout_horizontal)

        self.dialog_layer_arithmetic = LayerArithmeticDialog()

        # Add tool tray buttons
        self.button_layer_slice = self.add_tool_button(
            description='Create layer slice',
            icon_path=os.path.join(ICON_PATH, "Stanley Knife-48.png"),
            category='transformations',
            enabled=False,
            callback=lambda: self.add_layer(
                window=self.active_window, layer=self.current_layer,
                from_roi=True))

    def setup_connections(self):
        # -- Communications setup
        # Listen for layer selection events, enable/disable buttons
        self.tree_widget_layer_list.itemSelectionChanged.connect(
            lambda: self.toggle_buttons(self.current_layer_item))

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
        self.button_remove_layer.clicked.connect(lambda:
                                                 self.remove_layer_item())

        # When the arithmetic button is clicked, show math dialog
        self.button_layer_arithmetic.clicked.connect(
            self._show_arithmetic_dialog)

        # Create a new layer based on any active ROIs
        # self.button_create_layer_slice.clicked.connect(
        #     lambda: Dispatch.on_add_roi_layer.emit(layer=self.current_layer,
        #                                            from_roi=True))

        # Allow changing of plot color
        self.button_change_color.clicked.connect(
            self._change_plot_color
        )

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

    @DispatchHandle.register_listener("on_add_layer")
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

    @DispatchHandle.register_listener("on_remove_layer")
    def remove_layer_item(self, layer=None):
        if layer is None:
            layer = self.current_layer

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

        Dispatch.on_removed_layer.emit(layer=layer, window=self.active_window)

    @DispatchHandle.register_listener("on_added_plot", "on_updated_plot")
    def update_layer_item(self, plot=None, *args, **kwargs):
        if plot is None:
            return

        layer = plot._layer
        pixmap = QPixmap(10, 10)
        pixmap.fill(plot.pen.color())
        icon = QIcon(pixmap)

        layer_item = self.get_layer_item(layer)

        if layer_item is not None:
            layer_item.setIcon(0, icon)

    def add_layer(self, layer=None, mask=None, window=None, from_roi=True):
        """
        Creates a layer object from the current ROIs of the active plot layer.

        Parameters
        ----------
        layer : specviz.core.data.Layer
            The current active layer of the active plot.
        window : QtGui.QMdiSubWindow
            The parent object within which the plot window resides.
        mask : ndarray
            Boolean mask.
        """
        # User attempts to slice before opening a file
        if layer is None and window is None:
            logging.error(
                "Cannot add new layer; no layer and no window "
                "provided.")
            return

        roi_mask = mask if mask is not None and not from_roi else \
            window.get_roi_mask(layer=layer)

        new_layer = Layer(layer._source, mask=roi_mask,
                          name=layer._source.name + "Layer Slice")

        Dispatch.on_add_layer.emit(layer=new_layer)

    def _show_arithmetic_dialog(self):
        if self.current_layer is None:
            return

        if self.dialog_layer_arithmetic.exec_():
            formula = self.dialog_layer_arithmetic\
                .line_edit_formula.text()

            current_window = self.active_window
            current_layers = self.all_layers
            new_layer = Layer.from_formula(formula, layers=current_layers)

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
                Dispatch.on_add_layer.emit(layer=new_layer)
            else:
                logging.info("{} not equivalent to {}.".format(
                    new_layer.data.unit, current_window._plot_units[1]))
                Dispatch.on_add_window(layer=new_layer)

    def _change_plot_color(self):
        plot = self.active_window.get_plot(self.current_layer)

        col = QColorDialog.getColor(
            plot._pen_stash['pen_on'].color(),
            self.tree_widget_layer_list)

        if col.isValid():
            plot.pen = col

            Dispatch.on_updated_plot.emit(container=plot)
        else:
            logging.warning("Color is not valid.")

    def toggle_buttons(self, layer_item):
        if layer_item is not None:
            self.button_layer_arithmetic.setEnabled(True)
            self.button_remove_layer.setEnabled(True)
            self.button_layer_slice.setEnabled(True)
            self.button_change_color.setEnabled(True)
        else:
            self.button_layer_arithmetic.setEnabled(False)
            self.button_remove_layer.setEnabled(False)
            self.button_layer_slice.setEnabled(False)
            self.button_change_color.setEnabled(False)

    @DispatchHandle.register_listener("on_activated_window")
    def toggle_enabled(self, window):
        if window is None:
            return

        layers = window.get_all_layers()

        self.tree_widget_layer_list.clear()

        for layer in layers:
            self.add_layer_item(layer)
            plot = window.get_plot(layer)
            self.update_layer_item(plot)

    @DispatchHandle.register_listener("on_clicked_layer")
    def _set_layer_visibility(self, layer_item, col=0):
        """
        Toggles the visibility of the plot in the sub window.

        Parameters
        ----------
        layer : Layer
            Layer object to toggle visibility.

        col : int
            QtTreeWidget data column.
        """
        layer = layer_item.data(0, Qt.UserRole)

        if layer is not None:
            current_window = self.active_window

            current_window.set_visibility(
                layer, layer_item.checkState(col) == Qt.Checked,
                override=True)


