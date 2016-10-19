from ...third_party.qtpy.QtWidgets import *
from ...third_party.qtpy.QtCore import *
from ...third_party.qtpy.QtGui import *
from ...core.comms import dispatch, DispatchHandle
from ...core.data import Spectrum1DRefLayer
from .sub_windows import PlotSubWindow


class UiMainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(UiMainWindow, self).__init__(parent)
        DispatchHandle.setup(self)

        self.showMaximized()
        self.setMinimumSize(QSize(640, 480))
        self.setDockOptions(QMainWindow.AnimatedDocks)
        self.setWindowTitle("SpecViz")

        self.widget_central = QWidget(self)
        self.setCentralWidget(self.widget_central)

        # Toolbar
        self.layout_vertical = QVBoxLayout(self.widget_central)

        # MDI area setup
        self.mdi_area = MdiArea(self.widget_central)
        self.mdi_area.setFrameShape(QFrame.StyledPanel)
        self.mdi_area.setFrameShadow(QFrame.Plain)
        self.mdi_area.setLineWidth(2)
        brush = QBrush(QColor(200, 200, 200))
        brush.setStyle(Qt.SolidPattern)
        self.mdi_area.setBackground(brush)
        self.mdi_area.setAcceptDrops(True)

        self.layout_vertical.addWidget(self.mdi_area)

        # Menu bar setup
        self.menu_bar = QMenuBar(self)

        self.menu_file = QMenu(self.menu_bar)
        self.menu_file.setTitle("File")
        self.menu_edit = QMenu(self.menu_bar)
        self.menu_edit.setTitle("Edit")
        self.menu_view = QMenu(self.menu_bar)
        self.menu_edit.setTitle("View")

        self.menu_docks = QMenu(self.menu_bar)

        self.setMenuBar(self.menu_bar)

        # Status bar setup
        self.status_bar = QStatusBar(self)

        self.setStatusBar(self.status_bar)


class MainWindow(UiMainWindow):
    def __init__(self, parent=None, *args, **kwargs):
        super(MainWindow, self).__init__(parent)

        self.mdi_area.subWindowActivated.connect(self._set_activated_window)

    def _set_activated_window(self, window):
        if window is None:
            all_sws = self.mdi_area.subWindowList(order=self.mdi_area.ActivationHistoryOrder)

            if len(all_sws) > 0:
                window = all_sws[-1]
            else:
                window = None

        dispatch.on_activated_window.emit(
            window=window.widget() if window is not None else None)

    @DispatchHandle.register_listener("on_add_window")
    def add_sub_window(self, data=None, layer=None, window=None, *args,
                       **kwargs):
        layer = layer or Spectrum1DRefLayer.from_parent(data)
        window = window or PlotSubWindow()
        window.add_plot(layer=layer)

        if window is not None:
            mdi_sub_window = self.mdi_area.addSubWindow(window)
            window.show()
            self._set_activated_window(mdi_sub_window)

    @DispatchHandle.register_listener("on_add_roi")
    def add_roi(self, *args, **kwargs):
        mdi_sub_window = self.mdi_area.activeSubWindow()
        window = mdi_sub_window.widget()
        window.add_roi()

    @DispatchHandle.register_listener("on_status_message")
    def update_message(self, message, timeout=0):
        self.status_bar.showMessage(message, timeout)


class MdiArea(QMdiArea):
    def __init__(self, *args, **kwargs):
        super(MdiArea, self).__init__(*args, **kwargs)

    def dragEnterEvent(self, e):
        if True:
            e.accept()
        else:
            e.ignore()

    def dropEvent(self, e):
        dispatch.on_add_window.emit(data=e.mimeData.data())
