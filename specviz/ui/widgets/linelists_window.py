from ...third_party.qtpy.QtWidgets import *
from ...third_party.qtpy.QtGui import *
from ...third_party.qtpy.QtCore import *

from ...core.comms import Dispatch, DispatchHandle


#TODO work in progress

# The line list window must be a full fledged window and not a dialog.
# Dialogs do not support things like menu bars and central widgets.
# They are also a bit cumbersome to use when modal behavior is of no
# importance. Lets try to treat this as a window for now, and see how
# it goes.

class UiLinelistsWindow(object):

    # this code was taken as-is from the Designer.
    # Cleaning it up sounds like a lower priority
    # task for now.
    def setupUi(self, MainWindow, title):
        MainWindow.setWindowTitle(title)
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(500, 500)
        MainWindow.setMinimumSize(QSize(300, 350))
        self.centralWidget = QWidget(MainWindow)
        self.centralWidget.setObjectName("centralWidget")
        self.gridLayout = QGridLayout(self.centralWidget)
        self.gridLayout.setContentsMargins(11, 11, 11, 11)
        self.gridLayout.setSpacing(6)
        self.gridLayout.setObjectName("gridLayout")
        self.horizontalLayout_5 = QHBoxLayout()
        self.horizontalLayout_5.setContentsMargins(11, 11, 11, 11)
        self.horizontalLayout_5.setSpacing(6)
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.lines_selected_label = QLabel(self.centralWidget)
        self.lines_selected_label.setObjectName("lines_selected_label")
        self.horizontalLayout_5.addWidget(self.lines_selected_label)
        self.label = QLabel(self.centralWidget)
        self.label.setObjectName("label")
        self.horizontalLayout_5.addWidget(self.label)
        self.draw_button = QPushButton(self.centralWidget)
        self.draw_button.setObjectName("draw_button")
        self.horizontalLayout_5.addWidget(self.draw_button)
        self.erase_button = QPushButton(self.centralWidget)
        self.erase_button.setObjectName("erase_button")
        self.horizontalLayout_5.addWidget(self.erase_button)
        self.dismiss_button = QPushButton(self.centralWidget)
        self.dismiss_button.setObjectName("dismiss_button")
        self.horizontalLayout_5.addWidget(self.dismiss_button)
        self.gridLayout.addLayout(self.horizontalLayout_5, 4, 0, 1, 1)
        self.verticalLayout_11 = QVBoxLayout()
        self.verticalLayout_11.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout_11.setSpacing(6)
        self.verticalLayout_11.setObjectName("verticalLayout_11")
        self.tabWidget = QTabWidget(self.centralWidget)
        self.tabWidget.setObjectName("tabWidget")
        self.verticalLayout_11.addWidget(self.tabWidget)
        self.gridLayout.addLayout(self.verticalLayout_11, 0, 0, 1, 1)
        self.horizontalLayout_7 = QHBoxLayout()
        self.horizontalLayout_7.setContentsMargins(11, 11, 11, 11)
        self.horizontalLayout_7.setSpacing(6)
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.add_set_button = QPushButton(self.centralWidget)
        self.add_set_button.setObjectName("add_set_button")
        self.horizontalLayout_7.addWidget(self.add_set_button)
        spacerItem = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
        self.horizontalLayout_7.addItem(spacerItem)
        self.gridLayout.addLayout(self.horizontalLayout_7, 2, 0, 2, 1)
        MainWindow.setCentralWidget(self.centralWidget)
        self.menuBar = QMenuBar(MainWindow)
        self.menuBar.setGeometry(QRect(0, 0, 767, 22))
        self.menuBar.setObjectName("menuBar")
        self.menuFile = QMenu(self.menuBar)
        self.menuFile.setObjectName("menuFile")
        MainWindow.setMenuBar(self.menuBar)
        self.mainToolBar = QToolBar(MainWindow)
        self.mainToolBar.setMovable(False)
        self.mainToolBar.setFloatable(False)
        self.mainToolBar.setObjectName("mainToolBar")
        MainWindow.addToolBar(Qt.TopToolBarArea, self.mainToolBar)
        self.statusBar = QStatusBar(MainWindow)
        self.statusBar.setObjectName("statusBar")
        MainWindow.setStatusBar(self.statusBar)
        self.actionOpen = QAction(MainWindow)
        icon = QIcon()
        icon.addPixmap(QPixmap(":/img/Open Folder-48.png"), QIcon.Normal, QIcon.Off)
        self.actionOpen.setIcon(icon)
        self.actionOpen.setObjectName("actionOpen")
        self.actionExit = QAction(MainWindow)
        self.actionExit.setObjectName("actionExit")
        self.actionRemove = QAction(MainWindow)
        self.actionRemove.setObjectName("actionRemove")
        self.actionChange_Color = QAction(MainWindow)
        self.actionChange_Color.setObjectName("actionChange_Color")
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionExit)
        self.menuBar.addAction(self.menuFile.menuAction())
        self.mainToolBar.addAction(self.actionOpen)
        self.mainToolBar.addSeparator()

        self.retranslateUi(MainWindow)
        QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QCoreApplication.translate
        self.lines_selected_label.setText(_translate("MainWindow", "0"))
        self.label.setText(_translate("MainWindow", "lines selected"))
        self.draw_button.setText(_translate("MainWindow", "Draw"))
        self.erase_button.setText(_translate("MainWindow", "Erase"))
        self.dismiss_button.setText(_translate("MainWindow", "Dismiss"))
        self.add_set_button.setText(_translate("MainWindow", "Add set"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.actionOpen.setText(_translate("MainWindow", "Open"))
        self.actionExit.setText(_translate("MainWindow", "Exit"))
        self.actionRemove.setText(_translate("MainWindow", "Remove"))
        self.actionRemove.setToolTip(_translate("MainWindow", "Removes the selected layer"))
        self.actionChange_Color.setText(_translate("MainWindow", "Change Color"))
        self.actionChange_Color.setToolTip(_translate("MainWindow", "Change the line color selected layer"))


class LineListsWindow(UiLinelistsWindow):
    def __init__(self, plot_window, parent=None):
        super(LineListsWindow, self).__init__()

        # Builds GUI
        self._main_window = QMainWindow()
        self.setupUi(self._main_window, str(plot_window))

        # Request that line lists be read from wherever are they sources.
        Dispatch.on_request_linelists.emit()

        self.buildViews(plot_window)

        # Connect buttons to appropriate signals.
        #
        # Note that, for the Draw operation, we have to pass the table views to
        # the handler, even though it would be better to handle the row selections
        # all in here for the sake of encapsulation. This is so because this class
        # is not a QWidget or one of its subclasses, thus it cannot implement a
        # DispatchHandle signal handler.
        self.draw_button.clicked.connect(lambda:Dispatch.on_plot_linelists.emit(table_views=self._table_views))
        self.erase_button.clicked.connect(Dispatch.on_erase_linelabels.emit)
        self.dismiss_button.clicked.connect(Dispatch.on_dismiss_linelists_window.emit)

    def buildViews(self, plot_window):

        # Table views must be preserved in the instance so they can be
        # passed to whoever is going to do the actual line list plotting.
        # The plotting code must know which lines (table rows) are selected
        # in each line list.
        self._table_views = []

        for linelist in plot_window.linelists:

            table_model = LineListTableModel(linelist)

            if table_model.rowCount() > 0:
                table_view = QTableView()
                table_view.setModel(table_model)

                table_view.setSelectionBehavior(QAbstractItemView.SelectRows)
                table_view.horizontalHeader().setStretchLastSection(True)
                table_view.resizeColumnsToContents()
                comments = linelist.meta['comments']

                pane = self._buildLinelistPane(table_view, comments)

                self.tabWidget.addTab(pane, table_model.getName())

                self._table_views.append(table_view)

    def show(self):
        self._main_window.show()

    def hide(self):
        self._main_window.hide()

    def _buildLinelistPane(self, table, comments):
        pane = QWidget()

        layout = QVBoxLayout()
        layout.setSizeConstraint(QLayout.SetMaximumSize)
        info = QTextBrowser()
        info.setMaximumHeight(100)
        info.setAutoFillBackground(True)
        info.setStyleSheet("background-color: rgb(230,230,230);")

        for comment in comments:
            info.append(comment)

        layout.addWidget(info)
        layout.addWidget(table)
        pane.setLayout(layout)

        return pane


class LineListTableModel(QAbstractTableModel):

    def __init__(self, table, parent=None, *args):

        QAbstractTableModel.__init__(self, parent, *args)

        self._table = table

    def rowCount(self, index_parent=None, *args, **kwargs):
        return len(self._table.columns[0])

    def columnCount(self, index_parent=None, *args, **kwargs):
        return len(self._table.columns)

    def data(self, index, role=None):
        if not index.isValid():
            return QVariant()
        elif role != Qt.DisplayRole:
            return QVariant()
        return QVariant(str(self._table.columns[index.column()][index.row()]))

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return self._table.colnames[section]
        return QAbstractTableModel.headerData(self, section, orientation, role)

    def getName(self):
        return self._table.name
