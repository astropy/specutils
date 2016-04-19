from ...third_party.qtpy.QtWidgets import *
from ...third_party.qtpy.QtGui import *
from ...third_party.qtpy.QtCore import *

from ..qt.layer_arithmetic_dialog import Ui_LayerArithmeticDialog
from ..qt.unit_change_dialog import Ui_UnitChangeDialog


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
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel
                                          | QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.layout_vertical.addWidget(self.buttonBox)

        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)


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


class LayerArithmeticDialog(QDialog):
    def __init__(self, parent=None):
        super(LayerArithmeticDialog, self).__init__(parent)
        # Run the widget setup
        self.ui_layer_arithmetic_dialog = Ui_LayerArithmeticDialog()
        self.ui_layer_arithmetic_dialog.setupUi(self)


class UnitChangeDialog(QDialog):
    def __init__(self, parent=None):
        super(UnitChangeDialog, self).__init__(parent)

        self.flux_unit = ''
        self.disp_unit = ''

        # Run the widget setup
        self.ui_unit_change_dialog = Ui_UnitChangeDialog()
        self.ui_unit_change_dialog.setupUi(self)

    def accept(self):
        self.flux_unit = self.ui_unit_change_dialog.fluxUnitLineEdit.text()
        self.disp_unit = self.ui_unit_change_dialog.dispersionUnitLineEdit.text()

        super(UnitChangeDialog, self).accept()

    def reject(self):
        super(UnitChangeDialog, self).reject()
