from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..third_party.qtpy.QtCore import *
from ..third_party.qtpy.QtWidgets import *
from ..third_party.qtpy.QtGui import *

from .qt.mainwindow import Ui_MainWindow
from .widgets.sub_windows import PlotSubWindow
from .widgets.dialogs import LayerArithmeticDialog
from .widgets.menus import LayerContextMenu
from ..core.comms import Dispatch, DispatchHandle
from ..core.annotation import LineIDMarker


class Viewer(QMainWindow):
    """
    The `Viewer` is the main construction area for all GUI widgets. This
    object does **not** control the interactions between the widgets,
    but only their creation and placement.
    """
    def __init__(self, parent=None):
        super(Viewer, self).__init__(parent)
        self._current_sub_window = None

        self.main_window = Ui_MainWindow()
        self.main_window.setupUi(self)
        self.wgt_data_list = self.main_window.listWidget
        self.wgt_layer_list = self.main_window.treeWidget_2
        self.wgt_model_list = self.main_window.treeWidget
        self.wgt_model_list.setHeaderLabels(["Parameter", "Value"])

        # Setup
        self._setup_connections()

        # Context menus
        self.wgt_layer_list.setContextMenuPolicy(Qt.CustomContextMenu)
        self.layer_context_menu = LayerContextMenu()

        self.wgt_model_list.setContextMenuPolicy(Qt.CustomContextMenu)

        # Define the layer arithmetic dialog
        self._layer_arithmetic_dialog = LayerArithmeticDialog()

        # Setup event handler
        DispatchHandle.setup(self)

    def _setup_connections(self):
        # Listen for subwindow selection events, update layer list on selection
        self.main_window.mdiArea.subWindowActivated.connect(
            self._set_current_sub_window)

        # When a user edits the model parameter field, validate the input
        self.wgt_model_list.itemChanged.connect(
                self._model_parameter_validation)

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

    def _set_current_sub_window(self, sub_window):
        sub_window = sub_window or self.main_window.mdiArea.currentSubWindow()

        if sub_window is None:
            sub_window = self.main_window.mdiArea.activatePreviousSubWindow()

        if self._current_sub_window != sub_window:
            self._current_sub_window = sub_window
            Dispatch.on_selected_window.emit(window=self._current_sub_window)

    @property
    def current_data(self):
        """
        Returns the currently selected data object from the data list widget.

        Returns
        -------
        data : specviz.core.data.Data
            The `Data` object of the currently selected row.
        """
        data_item = self.wgt_data_list.currentItem()

        if data_item is not None:
            data = data_item.data(Qt.UserRole)
            return data

    @property
    def current_layer(self):
        """
        Returns the currently selected layer object form the layer list widget.

        Returns
        -------
        layer : specviz.core.data.Layer
            The `Layer` object of the currently selected row.
        """
        layer_item = self.wgt_layer_list.currentItem()

        if layer_item is not None:
            layer = layer_item.data(0, Qt.UserRole)

            return layer

    @property
    def current_layer_item(self):
        return self.wgt_layer_list.currentItem()

    @property
    def current_model_item(self):
        return self.wgt_model_list.currentItem()

    @property
    def current_sub_window(self):
        """
        Returns the currently active `PlotSubWindow` object.

        Returns
        -------
        sub_window : PlotSubWindow
            The currently active `PlotSubWindow` object.
        """
        if self._current_sub_window is not None:
            return self._current_sub_window.widget()

    @property
    def current_model(self):
        return self.main_window.modelsComboBox.currentText()

    @property
    def current_fitter(self):
        return self.main_window.fittingRoutinesComboBox.currentText()

    @property
    def current_model_formula(self):
        return self.main_window.lineEdit.text()

    def add_sub_window(self, *args, **kwargs):
        """
        Creates a new sub window instance in the MDI area.

        Returns
        -------
        new_sub_window : QMdiSubWindow
            The MdiSubWindow Qt instance.
        wgt_sub_window : QWidget
            The widget object within the QMdiSubWindow.
        """
        # Create new window
        plot_sub_window = PlotSubWindow()

        new_sub_window = self.main_window.mdiArea.addSubWindow(plot_sub_window)
        new_sub_window.show()

        return plot_sub_window

    def open_file_dialog(self, filters):
        """
        Given a list of filters, prompts the user to select an existing file
        and returns the file path and filter.

        Parameters
        ----------
        filters : list
            List of filters for the dialog.

        Returns
        -------
        file_name : str
            Path to the selected file.
        selected_filter : str
            The chosen filter (this indicates which custom loader from the
            registry to use).
        """
        dialog = QFileDialog(self)
        dialog.setFileMode(QFileDialog.ExistingFile)
        dialog.setNameFilters([x for x in filters])

        if dialog.exec_():
            file_names = dialog.selectedFiles()
            selected_filter = dialog.selectedNameFilter()

            return file_names[0], selected_filter

        return None, None

    @DispatchHandle.register_listener("on_added_data")
    def add_data_item(self, data):
        """
        Adds a `Data` object to the loaded data list widget.

        Parameters
        ----------
        data : specviz.core.data.Data
            The `Data` object to add to the list widget.
        """
        new_item = QListWidgetItem(data.name, self.wgt_data_list)
        new_item.setFlags(new_item.flags() |  Qt.ItemIsEditable)

        new_item.setData(Qt.UserRole, data)

        self.wgt_data_list.setCurrentItem(new_item)

    @DispatchHandle.register_listener("on_removed_data")
    def remove_data_item(self, data):
        data_item = self.get_data_item(data)

        self.wgt_data_list.takeItem(self.wgt_data_list.row(data_item))

    def get_data_item(self, data):
        for i in range(self.wgt_data_list.count()):
            data_item = self.wgt_data_list.item(0)

            if data_item.data(Qt.UserRole) == data:
                return data_item

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

        new_item = QTreeWidgetItem(self.get_layer_item(layer._parent) or
                                   self.wgt_layer_list)
        new_item.setFlags(new_item.flags() | Qt.ItemIsUserCheckable | Qt.ItemIsEditable)
        new_item.setText(0, layer.name)
        new_item.setData(0, Qt.UserRole, layer)
        new_item.setCheckState(0, Qt.Checked)

        self.wgt_layer_list.setCurrentItem(new_item)

    def get_layer_item(self, layer):
        root = self.wgt_layer_list.invisibleRootItem()

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
        root = self.wgt_layer_list.invisibleRootItem()

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

    def clear_layer_widget(self):
        self.wgt_layer_list.clear()

    def clear_model_widget(self):
        self.wgt_model_list.clear()

    @DispatchHandle.register_listener("on_updated_stats")
    def update_statistics(self, stats, layer):
        self.main_window.currentLayerLineEdit.setText("{}".format(layer.name))

        if 'mean' in stats:
            self.main_window.meanLineEdit.setText("{0:4.4g}".format(
                stats['mean'].value))

            self.main_window.medianLineEdit.setText("{0:4.4g}".format(
                stats['median'].value))

            self.main_window.standardDeviationLineEdit.setText("{0:4.4g}".format(
                stats['stddev'].value))

            self.main_window.totalLineEdit.setText("{0:4.4g}".format(
                float(stats['total'].value)))

            self.main_window.dataPointCountLineEdit.setText(
                str(stats['npoints']))

        if 'eq_width' in stats:
            self.main_window.equivalentWidthLineEdit.setText("{0:4.4g}".format(
                float(stats['eq_width'].value)))

        if 'centroid' in stats:
            self.main_window.centroidLineEdit.setText("{0:5.5g}".format(
                float(stats['centroid'].value)))

        if 'flux' in stats:
            self.main_window.fluxLineEdit.setText("{0:4.4g}".format(
                float(stats['flux'].value)))

        if 'avg_cont' in stats:
            self.main_window.meanContinuumLineEdit.setText("{0:4.4g}".format(
                float(stats['avg_cont'].value)))


    @DispatchHandle.register_listener("on_added_linelist")
    def add_linelist(self, linelist):

        # This is setting all markers at a fixed heigth in the
        # initial (before any zoom) data coordinates. Still TBD
        # how to do this in the generic case. Maybe derive heights
        # from curve data instead? Make the markers follow the
        # curve ups and downs?
        #
        # Ideally we would like to have the marker's X coordinate
        # pinned down to the plot surface in data value, and the Y
        # coordinate pinned down in screen value. This would make
        # the markers to stay at the same height in the window even
        # when the plot is zoomed. This kind of functionality doesn't
        # seem to be possible under pyqtgraph though. This requires
        # more investigation.

        plot_item = self.current_sub_window._plot_item

        # curve = plot_item.curves[0]

        data_range = plot_item.vb.viewRange()
        ymin = data_range[1][0]
        ymax = data_range[1][1]
        height = (ymax - ymin) * 0.75 + ymin

        # column names are defined in the YAML files.
        wave_column = linelist.columns['wavelength']
        id_column = linelist.columns['id']

        for i in range(len(wave_column)):
            marker = LineIDMarker(id_column[i], plot_item, orientation='vertical')

            marker.setPos(wave_column[i], height)

            plot_item.addItem(marker)
            # plot_item.addItem(marker.arrow)

            plot_item.update()

        plot_item.update()



