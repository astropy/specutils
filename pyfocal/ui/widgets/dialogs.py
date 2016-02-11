from qtpy.QtWidgets import *
from qtpy.QtGui import *


class TopAxisDialog(QDialog):
    def __init__(self, parent=None):
        super(TopAxisDialog, self).__init__(parent)
        self.ref_wave = 0.0
        self.redshift = 0.0

        self.vb_layout_main = QVBoxLayout()
        self.setLayout(self.vb_layout_main)

        self._container_list = []

        self.wgt_display_axis = QComboBox()
        self.wgt_display_axis.addItems(["Redshifted Wavelength", "Velocity",
                                        "Channel"])
        self.wgt_display_axis.currentIndexChanged.connect(self._on_select)

        frm_select = QFormLayout()
        frm_select.addRow("Display axis:", self.wgt_display_axis)

        # Redshift parameters
        self.grp_redshift = QGroupBox("Redshift Parameters")
        self.wgt_redshift = QLineEdit()
        self.wgt_redshift.setValidator(QDoubleValidator())
        frm_redshift = QFormLayout()
        self.grp_redshift.setLayout(frm_redshift)
        frm_redshift.addRow("Amount:", self.wgt_redshift)
        self._container_list.append(self.grp_redshift)

        # Velocity parameters
        self.grp_vel = QGroupBox("Velocity Parameters")
        self.wgt_ref_wave_unit = QLabel("")
        self.wgt_ref_wave = QLineEdit()
        hb_ref_wave = QHBoxLayout()
        hb_ref_wave.addWidget(self.wgt_ref_wave)
        hb_ref_wave.addWidget(self.wgt_ref_wave_unit)
        self.wgt_ref_wave.setValidator(QDoubleValidator())
        frm_vel = QFormLayout()
        self.grp_vel.setLayout(frm_vel)
        frm_vel.addRow("Reference Wavelength:", hb_ref_wave)
        self._container_list.append(self.grp_vel)

        button_box = QDialogButtonBox(QDialogButtonBox.Ok |
                                      QDialogButtonBox.Cancel)
        button_box.accepted.connect(self._on_accept)
        button_box.rejected.connect(self._on_reject)

        self.vb_layout_main.addLayout(frm_select)
        self.vb_layout_main.addWidget(self.grp_redshift)
        self.vb_layout_main.addWidget(self.grp_vel)
        self.vb_layout_main.addWidget(button_box)

        self._on_select(0)

    def set_current_unit(self, unit):
        self.wgt_ref_wave_unit.setText(unit)

    def _on_select(self, index):
        for cntr in self._container_list:
            cntr.hide()

        if index < len(self._container_list):
            self._container_list[index].show()

    def _on_accept(self):
        self.mode = self.wgt_display_axis.currentIndex()

        rw_val = str(self.wgt_ref_wave.text())
        self.ref_wave = float(rw_val) if rw_val != '' else self.ref_wave
        rs = str(self.wgt_redshift.text())
        self.redshift = float(rs) if rs != '' else self.redshift

        super(TopAxisDialog, self).accept()

    def _on_reject(self):
        super(TopAxisDialog, self).reject()