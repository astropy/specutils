from ...third_party.qtpy.QtWidgets import *
from ...third_party.qtpy.QtGui import *

from ..qt.axisdialog import Ui_Dialog


class TopAxisDialog(QDialog):
    def __init__(self, parent=None):
        super(TopAxisDialog, self).__init__(parent)
        self.ref_wave = 0.0
        self.redshift = 0.0

        # Run the widget setup
        self.ui_axis_dialog = Ui_Dialog()
        self.ui_axis_dialog.setupUi(self)

        # Populate options
        self.ui_axis_dialog.axisModeComboBox.addItems(['Velocity',
                                                       'Redshift', 'Channel'])

        # Setup connections
        self._setup_connections()

    def _setup_connections(self):
        # Show/hide corresponding container when mode is selected
        self.ui_axis_dialog.axisModeComboBox.currentIndexChanged.connect(self._on_select)

    def set_current_unit(self, unit):
        self.wgt_ref_wave_unit.setText(unit)

    def _on_select(self, index):
        pass

    def accept(self):
        self.mode = self.ui_axis_dialog.axisModeComboBox.currentIndex()

        rw_val = str(self.wgt_ref_wave.text())
        self.ref_wave = float(rw_val) if rw_val != '' else self.ref_wave
        rs = str(self.wgt_redshift.text())
        self.redshift = float(rs) if rs != '' else self.redshift

        super(TopAxisDialog, self).accept()

    def reject(self):
        super(TopAxisDialog, self).reject()