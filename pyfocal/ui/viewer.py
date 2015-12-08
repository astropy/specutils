from __future__ import absolute_import, division, print_function

from qtpy.QtGui import *
from .qt.mainwindow import Ui_MainWindow
from .qt.spectrasubwindow import Ui_SpectraSubWindow


class Viewer(QMainWindow):
    def __init__(self, parent=None):
        super(Viewer, self).__init__(parent)
        self.main_window = Ui_MainWindow()
        self.main_window.setupUi(self)

        self.add_sub_window()

    def add_sub_window(self):
        wgt_sub_window = QMainWindow()
        ui_spectra_sub_window = Ui_SpectraSubWindow()
        ui_spectra_sub_window.setupUi(wgt_sub_window)
        new_sub_window = self.main_window.mdiArea.addSubWindow(wgt_sub_window)
        return new_sub_window, ui_spectra_sub_window.webView