from qtpy.QtGui import *
from qtpy.QtCore import *
from qtpy.QtWidgets import *
# from PyQt5.QtWidgets import *
from .qt.mainwindow import Ui_MainWindow
from .qt.plotsubwindow import Ui_SpectraSubWindow
from .widgets.plot_window import PlotWindow


class Viewer(QMainWindow):
    """
    The `Viewer` is the main construction area for all GUI widgets. This
    object does **not** control the interactions between the widgets,
    but only their creation and placement.
    """
    def __init__(self, parent=None):
        super(Viewer, self).__init__(parent)
        self.main_window = Ui_MainWindow()
        self.main_window.setupUi(self)
        self.wgt_data_list = self.main_window.listWidget
        self.wgt_layer_list = self.main_window.treeWidget_2
        self.wgt_model_list = self.main_window.treeWidget
        self.wgt_model_list.setHeaderLabels(["Parameter", "Value"])

        # Connect the validation events
        self.wgt_model_list.itemChanged.connect(
                self._model_parameter_validation)

        # Setup context menus
        self._setup_context_menus()

    def _setup_context_menus(self):
        self.wgt_layer_list.customContextMenuRequested.connect(
                self._layer_context_menu)

    def _set_model_tool_options(self):
        layer = self.current_layer()

        if layer is None:
            return

        if not hasattr(layer, 'model'):
            self.main_window.pushButton_4.show()
            self.main_window.pushButton_2.hide()
        else:
            self.main_window.pushButton_4.hide()
            self.main_window.pushButton_2.show()

    @property
    def current_model(self):
        return self.main_window.comboBox.currentText()

    @property
    def current_model_formula(self):
        return self.main_window.lineEdit.text()

    def add_sub_window(self):
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
        plot_sub_window = PlotWindow()

        # Populate window with tool bars, status, etc.
        ui_sub_window = Ui_SpectraSubWindow()
        ui_sub_window.setupUi(plot_sub_window)

        # Let the sub window do initialization
        plot_sub_window.initialize()

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

    def add_data_item(self, data):
        """
        Adds a `Data` object to the loaded data list widget.

        Parameters
        ----------
        data : pyfocal.core.data.Data
            The `Data` object to add to the list widget.
        """
        new_item = QListWidgetItem(data.name, self.wgt_data_list)
        new_item.setData(Qt.UserRole, data)

    def add_layer_item(self, layer, *args):
        """
        Adds a `Layer` object to the loaded layer list widget.

        Parameters
        ----------
        layer : pyfocal.core.data.Layer
            The `Layer` object to add to the list widget.
        """
        new_item = QTreeWidgetItem(self.get_layer_item(layer._source) or
                                   self.wgt_layer_list)
        new_item.setText(0, layer.name)
        new_item.setData(0, Qt.UserRole, layer)

        self.wgt_layer_list.setCurrentItem(new_item)

    def get_layer_item(self, layer):
        root = self.wgt_layer_list.invisibleRootItem()

        for i in range(root.childCount()):
            child = root.child(i)

            if child.data(0, Qt.UserRole) == layer:
                return child

    def remove_layer_item(self, layer):
        for child in self.wgt_layer_list.children():
            if child.data(Qt.UserRole) == layer:
                self.wgt_layer_list.removeItemWidget(child)
                break

    def add_model_item(self, model):
        """
        Adds an `astropy.modeling.Model` to the loaded model tree widget.

        Parameters
        ----------
        """
        name = model.__class__.__name__
        new_item = QTreeWidgetItem(self.wgt_model_list)
        new_item.setFlags(new_item.flags() | Qt.ItemIsEditable)

        new_item.setText(0, name)
        new_item.setData(0, Qt.UserRole, model)

        for i, para in enumerate(model.param_names):
            new_para_item = QTreeWidgetItem(new_item)
            new_para_item.setText(0, para)
            new_para_item.setData(0, Qt.UserRole,
                                  model.parameters[i])
            new_para_item.setText(1, str(model.parameters[i]))
            new_para_item.setFlags(new_para_item.flags() | Qt.ItemIsEditable)

    def remove_model_item(self, layer, model):
        for child in self.wgt_model_list.children():
            if child.data(Qt.UserRole) == model:
                self.wgt_model_list.removeItemWidget(child)
                break

    def _model_parameter_validation(self, item, col):
        if col == 0:
            return

        try:
            item.setText(col, str(float(item.text(col))))
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

    def current_data(self):
        """
        Returns the currently selected data object from the data list widget.

        Returns
        -------
        data : pyfocal.core.data.Data
            The `Data` object of the currently selected row.
        """
        data_item = self.wgt_data_list.currentItem()

        if data_item is not None:
            data = data_item.data(Qt.UserRole)
            return data

    def current_layer(self):
        """
        Returns the currently selected layer object form the layer list widget.

        Returns
        -------
        layer : pyfocal.core.data.Layer
            The `Layer` object of the currently selected row.
        """
        layer_item = self.wgt_layer_list.currentItem()

        if layer_item is not None:
            layer = layer_item.data(0, Qt.UserRole)

            return layer

    def current_sub_window(self):
        """
        Returns the currently active `QMdiSubWindow` object.

        Returns
        -------
        sub_window : QMdiSubWindow
            The currently active `QMdiSubWindow` object.
        """
        sub_window = self.main_window.mdiArea.currentSubWindow()

        if sub_window is not None:
            return sub_window.widget()

    def update_statistics(self, stat_dict):
        self.main_window.label_2.setText("{0:0.03f}".format(stat_dict['mean']))
        self.main_window.label_4.setText("{0:0.03f}".format(stat_dict['median']))
        self.main_window.label_6.setText("{0:0.03f}".format(stat_dict['stddev']))
        self.main_window.label_8.setText("{0:0.03f}".format(stat_dict['total']))
        self.main_window.label_10.setText(str(stat_dict['npoints']))

    def _layer_context_menu(self, point):
        menu = QMenu()
        menu.addAction(self.main_window.actionChange_Color)
        menu.addAction(self.main_window.actionRemove)
        menu.exec_(self.wgt_layer_list.viewport().mapToGlobal(point))
