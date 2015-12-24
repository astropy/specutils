from __future__ import absolute_import, division, print_function

from PyQt4.QtGui import *
from PyQt4.QtCore import *
# from PyQt5.QtWidgets import *
from .qt.mainwindow import Ui_MainWindow
from .qt.spectrasubwindow import Ui_SpectraSubWindow


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
        self.wgt_layer_list = self.main_window.listWidget_2

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
        main_window = QMainWindow()
        wgt_sub_window = Ui_SpectraSubWindow()
        wgt_sub_window.setupUi(main_window)
        new_sub_window = self.main_window.mdiArea.addSubWindow(main_window)

        return new_sub_window, wgt_sub_window

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
            selected_filter = dialog.selectedFilter()

            return file_names[0], selected_filter

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

    def add_layer_item(self, layer):
        """
        Adds a `Layer` object to the loaded layer list widget.

        Parameters
        ----------
        layer : pyfocal.core.data.Layer
            The `Layer` object to add to the list widget.
        """
        new_item = QListWidgetItem(layer.name, self.wgt_layer_list)
        new_item.setData(Qt.UserRole, layer)

    def current_data(self):
        """
        Returns the currently selected data object from the data list widget.

        Returns
        -------
        data : pyfocal.core.data.Data
            The `Data` object of the currently selected row.

        """
        data_item = self.wgt_data_list.currentItem()
        data = data_item.data(Qt.UserRole).toPyObject()

        return data
