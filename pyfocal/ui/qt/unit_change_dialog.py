# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file './source/unit_change_dialog.ui'
#
# Created by: ...third_party.qtpy UI code generator 5.5
#
# WARNING! All changes made in this file will be lost!
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...third_party.qtpy import QtCore, QtGui, QtWidgets

class Ui_UnitChangeDialog(object):
    def setupUi(self, UnitChangeDialog):
        UnitChangeDialog.setObjectName("UnitChangeDialog")
        UnitChangeDialog.resize(320, 118)
        self.verticalLayout = QtWidgets.QVBoxLayout(UnitChangeDialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setObjectName("formLayout")
        self.fluxUnitLabel = QtWidgets.QLabel(UnitChangeDialog)
        self.fluxUnitLabel.setObjectName("fluxUnitLabel")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.fluxUnitLabel)
        self.fluxUnitLineEdit = QtWidgets.QLineEdit(UnitChangeDialog)
        self.fluxUnitLineEdit.setObjectName("fluxUnitLineEdit")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.fluxUnitLineEdit)
        self.dispersionUnitLabel = QtWidgets.QLabel(UnitChangeDialog)
        self.dispersionUnitLabel.setObjectName("dispersionUnitLabel")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.dispersionUnitLabel)
        self.dispersionUnitLineEdit = QtWidgets.QLineEdit(UnitChangeDialog)
        self.dispersionUnitLineEdit.setObjectName("dispersionUnitLineEdit")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.dispersionUnitLineEdit)
        self.verticalLayout.addLayout(self.formLayout)
        self.buttonBox = QtWidgets.QDialogButtonBox(UnitChangeDialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(UnitChangeDialog)
        self.buttonBox.accepted.connect(UnitChangeDialog.accept)
        self.buttonBox.rejected.connect(UnitChangeDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(UnitChangeDialog)

    def retranslateUi(self, UnitChangeDialog):
        _translate = QtCore.QCoreApplication.translate
        UnitChangeDialog.setWindowTitle(_translate("UnitChangeDialog", "Change Plot Units"))
        self.fluxUnitLabel.setText(_translate("UnitChangeDialog", "Flux Unit"))
        self.dispersionUnitLabel.setText(_translate("UnitChangeDialog", "Dispersion Unit"))

