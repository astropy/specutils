from ...third_party.qtpy.QtWidgets import *
from ...third_party.qtpy.QtGui import *

from ..qt.axisdialog import Ui_Dialog
from ..qt.layer_arithmetic_dialog import Ui_LayerArithmeticDialog
from ..qt.unit_change_dialog import Ui_UnitChangeDialog


class TopAxisDialog(QDialog):
    def __init__(self, parent=None):
        super(TopAxisDialog, self).__init__(parent)
        self.ref_wave = 0.0
        self.redshift = 0.0

        # Run the widget setup
        self.ui_axis_dialog = Ui_Dialog()
        self.ui_axis_dialog.setupUi(self)

        # Set validators
        self.ui_axis_dialog.referenenceWavelengthLineEdit.setValidator(
            QDoubleValidator())
        self.ui_axis_dialog.redshiftAmountLineEdit.setValidator(
            QDoubleValidator())

        # Populate options
        self.ui_axis_dialog.axisModeComboBox.addItems(
            ['Velocity', 'Redshift', 'Pixel'])

        # Setup connections
        self._setup_connections()
        self._on_select(0)

    def _setup_connections(self):
        # Show/hide corresponding container when mode is selected
        self.ui_axis_dialog.axisModeComboBox.currentIndexChanged.connect(self._on_select)

    def set_current_unit(self, unit):
        self.wgt_ref_wave_unit.setText(unit)

    def _on_select(self, index):
        if index == 0:
            self.ui_axis_dialog.velocityGroupBox.show()
            self.ui_axis_dialog.redshiftGroupBox.hide()
        elif index == 1:
            self.ui_axis_dialog.velocityGroupBox.hide()
            self.ui_axis_dialog.redshiftGroupBox.show()
        else:
            self.ui_axis_dialog.velocityGroupBox.hide()
            self.ui_axis_dialog.redshiftGroupBox.hide()

    def accept(self):
        self.mode = self.ui_axis_dialog.axisModeComboBox.currentIndex()

        rw_val = str(self.ui_axis_dialog.referenenceWavelengthLineEdit.text())
        self.ref_wave = float(rw_val) if rw_val != '' else self.ref_wave
        rs = str(self.ui_axis_dialog.redshiftAmountLineEdit.text())
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
        self.disp_unit = self.ui_unit_change_dialog.dispersionUnitLineEdit\
            .text()

        super(UnitChangeDialog, self).accept()

    def reject(self):
        super(UnitChangeDialog, self).reject()
