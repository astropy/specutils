# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file './source/layer_arithmetic_dialog.ui'
#
# Created by: ...third_party.qtpy UI code generator 5.5.1
#
# WARNING! All changes made in this file will be lost!
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...third_party.qtpy import QtCore, QtGui, QtWidgets

class Ui_LayerArithmeticDialog(object):
    def setupUi(self, LayerArithmeticDialog):
        LayerArithmeticDialog.setObjectName("LayerArithmeticDialog")
        LayerArithmeticDialog.resize(354, 134)
        self.verticalLayout = QtWidgets.QVBoxLayout(LayerArithmeticDialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.layer1GroupBox = QtWidgets.QGroupBox(LayerArithmeticDialog)
        self.layer1GroupBox.setObjectName("layer1GroupBox")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.layer1GroupBox)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.formulaLineEdit = QtWidgets.QLineEdit(self.layer1GroupBox)
        self.formulaLineEdit.setObjectName("formulaLineEdit")
        self.verticalLayout_2.addWidget(self.formulaLineEdit)
        self.horizontalLayout.addWidget(self.layer1GroupBox)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.buttonBox = QtWidgets.QDialogButtonBox(LayerArithmeticDialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(LayerArithmeticDialog)
        self.buttonBox.accepted.connect(LayerArithmeticDialog.accept)
        self.buttonBox.rejected.connect(LayerArithmeticDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(LayerArithmeticDialog)

    def retranslateUi(self, LayerArithmeticDialog):
        _translate = QtCore.QCoreApplication.translate
        LayerArithmeticDialog.setWindowTitle(_translate("LayerArithmeticDialog", "Layer Arithmetic"))
        self.layer1GroupBox.setTitle(_translate("LayerArithmeticDialog", "Formula"))

