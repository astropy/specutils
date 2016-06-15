from ...third_party.qtpy.QtWidgets import *
from ...third_party.qtpy.QtGui import *
from ...third_party.qtpy.QtCore import *

from ...core.comms import Dispatch

class UiTopAxisDialog(QDialog):
    """
    Initialize all the TopAxisDialog Qt UI elements.
    """
    def __init__(self, *args, **kwargs):
        super(UiTopAxisDialog, self).__init__(*args, **kwargs)
        self.setObjectName("Top Axis Dialog")

        size_policy = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        size_policy.setHorizontalStretch(0)
        size_policy.setVerticalStretch(0)
        self.setSizePolicy(size_policy)

        # Dialog settings
        self.setWindowTitle("Axis Settings")

        self.layout_vertical = QVBoxLayout(self)
        self.layout_horizontal = QHBoxLayout()
        self.layout_vertical.addLayout(self.layout_horizontal)

        # Define header selectors
        self.label_axis_mode = QLabel(self)
        self.combo_box_axis_mode = QComboBox(self)

        self.label_axis_mode.setText("Axis mode")

        self.layout_horizontal.addWidget(self.label_axis_mode)
        self.layout_horizontal.addWidget(self.combo_box_axis_mode)

        # Define velocity
        self.group_box_velocity = QGroupBox(self)
        self.label_reference_wavelength = QLabel(self.group_box_velocity)
        self.line_edit_reference_wavelength = QLineEdit(self.group_box_velocity)

        self.group_box_velocity.setTitle("Velocity parameters")
        self.label_reference_wavelength.setText("Reference wavelength")

        self.layout_horizontal_2 = QHBoxLayout(self.group_box_velocity)
        self.layout_horizontal_2.addWidget(self.label_reference_wavelength)
        self.layout_horizontal_2.addWidget(self.line_edit_reference_wavelength)

        self.layout_vertical.addWidget(self.group_box_velocity)

        # Define redshift
        self.group_box_redshift = QGroupBox(self)
        self.label_redshift = QLabel(self.group_box_redshift)
        self.line_edit_redshift = QLineEdit(
            self.group_box_redshift)

        self.group_box_redshift.setTitle("Redshift parameters")
        self.label_redshift.setText("Amount")

        self.layout_horizontal_3 = QHBoxLayout(self.group_box_redshift)
        self.layout_horizontal_3.addWidget(self.label_redshift)
        self.layout_horizontal_3.addWidget(self.line_edit_redshift)

        self.layout_vertical.addWidget(self.group_box_redshift)

        # Add a spacer
        self.layout_vertical.addStretch(1)

        # Buttons
        self.button_box = QDialogButtonBox(self)
        self.button_box.setOrientation(Qt.Horizontal)
        self.button_box.setStandardButtons(QDialogButtonBox.Cancel
                                           | QDialogButtonBox.Ok)
        self.button_box.setObjectName("buttonBox")
        self.layout_vertical.addWidget(self.button_box)

        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)


class TopAxisDialog(UiTopAxisDialog):
    def __init__(self, parent=None):
        super(TopAxisDialog, self).__init__(parent)
        self.ref_wave = 0.0
        self.redshift = 0.0

        # Set validators
        self.line_edit_reference_wavelength.setValidator(
            QDoubleValidator())
        self.line_edit_redshift.setValidator(
            QDoubleValidator())

        # Populate options
        self.combo_box_axis_mode.addItems(
            ['Velocity', 'Redshift', 'Pixel'])

        # Setup connections
        self._setup_connections()
        self._on_select(0)

    def _setup_connections(self):
        # Show/hide corresponding container when mode is selected
        self.combo_box_axis_mode.currentIndexChanged.connect(self._on_select)

    def set_current_unit(self, unit):
        self.wgt_ref_wave_unit.setText(unit)

    def _on_select(self, index):
        if index == 0:
            self.group_box_velocity.show()
            self.group_box_redshift.hide()
        elif index == 1:
            self.group_box_velocity.hide()
            self.group_box_redshift.show()
        else:
            self.group_box_velocity.hide()
            self.group_box_redshift.hide()

    def accept(self):
        self.mode = self.combo_box_axis_mode.currentIndex()

        rw_val = str(self.line_edit_reference_wavelength.text())
        self.ref_wave = float(rw_val) if rw_val != '' else self.ref_wave
        rs = str(self.line_edit_redshift.text())
        self.redshift = float(rs) if rs != '' else self.redshift

        super(TopAxisDialog, self).accept()

    def reject(self):
        super(TopAxisDialog, self).reject()


class UiLayerArithmeticDialog(QDialog):
    def __init__(self, parent=None):
        super(UiLayerArithmeticDialog, self).__init__(parent)

        # Dialog settings
        self.setWindowTitle("Layer Arithmetic")
        self.resize(354, 134)

        self.layout_vertical = QVBoxLayout(self)

        # Arithmetic group box
        self.group_box_arithmetic = QGroupBox(self)
        self.group_box_arithmetic.setTitle("Formula")

        self.line_edit_formula = QLineEdit(self.group_box_arithmetic)

        self.layout_horizontal = QHBoxLayout(self.group_box_arithmetic)
        self.layout_horizontal.addWidget(self.line_edit_formula)
        self.layout_vertical.addWidget(self.group_box_arithmetic)

        # Buttons
        self.button_box = QDialogButtonBox(self)
        self.button_box.setOrientation(Qt.Horizontal)
        self.button_box.setStandardButtons(
            QDialogButtonBox.Cancel | QDialogButtonBox.Ok)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)

        self.layout_vertical.addWidget(self.button_box)


class LayerArithmeticDialog(UiLayerArithmeticDialog):
    def __init__(self, parent=None):
        super(LayerArithmeticDialog, self).__init__(parent)


class UiUnitChangeDialog(QDialog):
    def __init__(self, parent=None):
        super(UiUnitChangeDialog, self).__init__(parent)

        # Dialog settings
        self.setWindowTitle("Change Plot Units")

        self.layout_vertical = QVBoxLayout(self)
        self.form_layout = QFormLayout()
        self.layout_vertical.addLayout(self.form_layout)

        # Flux unit
        self.label_flux_unit = QLabel(self)
        self.line_edit_flux_unit = QLineEdit(self)

        self.label_flux_unit.setText("Flux Unit")

        self.form_layout.addRow(self.label_flux_unit, self.line_edit_flux_unit)

        # Dispersion unit
        self.label_disp_unit = QLabel(self)
        self.line_edit_disp_unit = QLineEdit(self)

        self.label_disp_unit.setText("Dispersion Unit")

        self.form_layout.addRow(self.label_disp_unit, self.line_edit_disp_unit)

        self.button_box = QDialogButtonBox(self)
        self.button_box.setOrientation(Qt.Horizontal)
        self.button_box.setStandardButtons(
            QDialogButtonBox.Cancel | QDialogButtonBox.Ok)
        self.layout_vertical.addWidget(self.button_box)

        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)


class UnitChangeDialog(UiUnitChangeDialog):
    def __init__(self, parent=None):
        super(UnitChangeDialog, self).__init__(parent)

        self.flux_unit = ''
        self.disp_unit = ''

    def accept(self):
        self.flux_unit = self.line_edit_flux_unit.text()
        self.disp_unit = self.line_edit_disp_unit.text()

        super(UnitChangeDialog, self).accept()

    def reject(self):
        super(UnitChangeDialog, self).reject()


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
    def setupUi(self, MainWindow):
        MainWindow.setWindowTitle("Line Lists")
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(767, 791)
        MainWindow.setMinimumSize(QSize(640, 480))
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
        self.tab = QWidget()
        self.tab.setObjectName("tab")
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QWidget()
        self.tab_2.setObjectName("tab_2")
        self.tabWidget.addTab(self.tab_2, "")
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
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "Tab 1"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "Tab 2"))
        self.add_set_button.setText(_translate("MainWindow", "Add set"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.actionOpen.setText(_translate("MainWindow", "Open"))
        self.actionExit.setText(_translate("MainWindow", "Exit"))
        self.actionRemove.setText(_translate("MainWindow", "Remove"))
        self.actionRemove.setToolTip(_translate("MainWindow", "Removes the selected layer"))
        self.actionChange_Color.setText(_translate("MainWindow", "Change Color"))
        self.actionChange_Color.setToolTip(_translate("MainWindow", "Change the line color selected layer"))


class LineListsWindow(UiLinelistsWindow):
    def __init__(self, parent=None):
        super(LineListsWindow, self).__init__()

        self.main_window = QMainWindow()
        self.setupUi(self.main_window)

        # Bare-bones implementation that has only the Draw and Erase buttons
        # active. The Draw action consists in requesting the input of pre-canned
        # line lists. The ingestion routine directly calls for the display of the
        # associated line labels. In a real-world implementation, the ingestion
        # should be handled separately from the drawing action, so as to give the
        # user the opportunity to interact with the line lists.
        self.draw_button.clicked.connect(Dispatch.on_request_linelist.emit)

        # The Erase action erases all line IDs that are currently plotted.
        self.erase_button.clicked.connect(Dispatch.on_erase_linelabels.emit)

    def show(self):
        self.main_window.show()



