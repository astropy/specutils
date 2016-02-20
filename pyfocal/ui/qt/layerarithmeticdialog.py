# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file './source/layerarithmeticdialog.ui'
#
# Created by: ...third_party.qtpy UI code generator 5.5
#
# WARNING! All changes made in this file will be lost!
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...third_party.qtpy import QtCore, QtGui, QtWidgets

class Ui_LayerArithmeticDialog(object):
    def setupUi(self, LayerArithmeticDialog):
        LayerArithmeticDialog.setObjectName("LayerArithmeticDialog")
        LayerArithmeticDialog.resize(442, 133)
        self.verticalLayout = QtWidgets.QVBoxLayout(LayerArithmeticDialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.layer1GroupBox = QtWidgets.QGroupBox(LayerArithmeticDialog)
        self.layer1GroupBox.setObjectName("layer1GroupBox")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.layer1GroupBox)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.layer1ComboBox = QtWidgets.QComboBox(self.layer1GroupBox)
        self.layer1ComboBox.setObjectName("layer1ComboBox")
        self.verticalLayout_2.addWidget(self.layer1ComboBox)
        self.horizontalLayout.addWidget(self.layer1GroupBox)
        self.operatorGroupBox = QtWidgets.QGroupBox(LayerArithmeticDialog)
        self.operatorGroupBox.setObjectName("operatorGroupBox")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.operatorGroupBox)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.operatorComboBox = QtWidgets.QComboBox(self.operatorGroupBox)
        self.operatorComboBox.setObjectName("operatorComboBox")
        self.verticalLayout_4.addWidget(self.operatorComboBox)
        self.horizontalLayout.addWidget(self.operatorGroupBox)
        self.layer2GroupBox = QtWidgets.QGroupBox(LayerArithmeticDialog)
        self.layer2GroupBox.setObjectName("layer2GroupBox")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.layer2GroupBox)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.layer2ComboBox = QtWidgets.QComboBox(self.layer2GroupBox)
        self.layer2ComboBox.setObjectName("layer2ComboBox")
        self.verticalLayout_3.addWidget(self.layer2ComboBox)
        self.horizontalLayout.addWidget(self.layer2GroupBox)
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
        self.layer1GroupBox.setTitle(_translate("LayerArithmeticDialog", "First Layer"))
        self.operatorGroupBox.setTitle(_translate("LayerArithmeticDialog", "Operator"))
        self.layer2GroupBox.setTitle(_translate("LayerArithmeticDialog", "Second Layer"))

