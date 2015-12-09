from __future__ import absolute_import, division, print_function

from PyQt4.QtGui import *
from PyQt4.QtCore import *
from .qt.mainwindow import Ui_MainWindow
from .qt.spectrasubwindow import Ui_SpectraSubWindow


class Viewer(QMainWindow):
    def __init__(self, parent=None):
        super(Viewer, self).__init__(parent)
        self.main_window = Ui_MainWindow()
        self.main_window.setupUi(self)
        self.wgt_data_list = self.main_window.listWidget

    def add_sub_window(self):
        wgt_sub_window = QMainWindow()
        ui_spectra_sub_window = Ui_SpectraSubWindow()
        ui_spectra_sub_window.setupUi(wgt_sub_window)
        new_sub_window = self.main_window.mdiArea.addSubWindow(wgt_sub_window)
        return new_sub_window

    def open_file_dialog(self, filters):
        dialog = QFileDialog(self)
        dialog.setFileMode(QFileDialog.ExistingFile)
        dialog.setNameFilters([x for x in filters])

        if dialog.exec_():
            file_names = dialog.selectedFiles()
            selected_filter = dialog.selectedFilter()
            return file_names[0], selected_filter

    def add_data_item(self, data):
        print("Adding data item")
        new_item = QListWidgetItem(data.name, self.wgt_data_list)
        new_item.setData(Qt.UserRole, data)

    def current_data(self):
        data_item = self.wgt_data_list.currentItem()
        return data_item.data(Qt.UserRole).toPyObject()
