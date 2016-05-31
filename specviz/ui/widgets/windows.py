from ...third_party.qtpy.QtWidgets import *
from ...third_party.qtpy.QtCore import *
from ...third_party.qtpy.QtGui import *
from ...core.comms import Dispatch, DispatchHandle


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

        # Tool bar setup
        # self.tool_bar_main = QToolBar(self)
        # self.tool_bar_main.setMovable(False)
        # self.tool_bar_main.setFloatable(False)
        #
        # self.action_open = QAction(self)
        # icon_open= QIcon()
        # icon_open.addPixmap(QPixmap(":/img/Open Folder-48.png"))
        # self.action_open.setIcon(icon_open)
        # self.tool_bar_main.addAction(self.action_open)
        #
        # self.addToolBar(Qt.TopToolBarArea, self.tool_bar_main)

        # Status bar setup
        self.status_bar = QStatusBar(self)

        self.setStatusBar(self.status_bar)


class MainWindow(UiMainWindow):
    def __init__(self, parent=None, *args, **kwargs):
        super(MainWindow, self).__init__(parent)

        self.mdi_area.subWindowActivated.connect(
            lambda x: Dispatch.on_activated_window.emit(
                window=x.widget() if x is not None else None))

    @DispatchHandle.register_listener("on_added_window")
    def add_sub_window(self, window=None, *args, **kwargs):
        if window is not None:
            mdi_sub_window = self.mdi_area.addSubWindow(window)
            window.show()

    @DispatchHandle.register_listener("on_add_roi")
    def add_roi(self, *args, **kwargs):
        mdi_sub_window = self.mdi_area.activeSubWindow()
        window = mdi_sub_window.widget()
        window.add_roi()


class MdiArea(QMdiArea):
    def __init__(self, *args, **kwargs):
        super(MdiArea, self).__init__(*args, **kwargs)

    def dragEnterEvent(self, e):
        if True:
            e.accept()
        else:
            e.ignore()

    def dropEvent(self, e):
        Dispatch.on_add_window.emit(data=e.mimeData.data())
