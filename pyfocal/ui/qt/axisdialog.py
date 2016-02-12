# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file './source/axisdialog.ui'
#
# Created by: ...third_party.qtpy UI code generator 5.5.1
#
# WARNING! All changes made in this file will be lost!
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...third_party.qtpy import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(300, 224)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Dialog.sizePolicy().hasHeightForWidth())
        Dialog.setSizePolicy(sizePolicy)
        self.verticalLayout = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.axisModeLabel = QtWidgets.QLabel(Dialog)
        self.axisModeLabel.setObjectName("axisModeLabel")
        self.horizontalLayout.addWidget(self.axisModeLabel)
        self.axisModeComboBox = QtWidgets.QComboBox(Dialog)
        self.axisModeComboBox.setObjectName("axisModeComboBox")
        self.horizontalLayout.addWidget(self.axisModeComboBox)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.velocityGroupBox = QtWidgets.QGroupBox(Dialog)
        self.velocityGroupBox.setObjectName("velocityGroupBox")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.velocityGroupBox)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.referenenceWavelengthLabel = QtWidgets.QLabel(self.velocityGroupBox)
        self.referenenceWavelengthLabel.setObjectName("referenenceWavelengthLabel")
        self.horizontalLayout_3.addWidget(self.referenenceWavelengthLabel)
        self.referenenceWavelengthLineEdit = QtWidgets.QLineEdit(self.velocityGroupBox)
        self.referenenceWavelengthLineEdit.setObjectName("referenenceWavelengthLineEdit")
        self.horizontalLayout_3.addWidget(self.referenenceWavelengthLineEdit)
        self.verticalLayout_2.addLayout(self.horizontalLayout_3)
        self.verticalLayout.addWidget(self.velocityGroupBox)
        self.redshiftGroupBox = QtWidgets.QGroupBox(Dialog)
        self.redshiftGroupBox.setObjectName("redshiftGroupBox")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.redshiftGroupBox)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.referenenceWavelengthLabel_2 = QtWidgets.QLabel(self.redshiftGroupBox)
        self.referenenceWavelengthLabel_2.setObjectName("referenenceWavelengthLabel_2")
        self.horizontalLayout_2.addWidget(self.referenenceWavelengthLabel_2)
        self.referenenceWavelengthLineEdit_2 = QtWidgets.QLineEdit(self.redshiftGroupBox)
        self.referenenceWavelengthLineEdit_2.setObjectName("referenenceWavelengthLineEdit_2")
        self.horizontalLayout_2.addWidget(self.referenenceWavelengthLineEdit_2)
        self.verticalLayout_3.addLayout(self.horizontalLayout_2)
        self.verticalLayout.addWidget(self.redshiftGroupBox)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(Dialog)
        self.buttonBox.accepted.connect(Dialog.accept)
        self.buttonBox.rejected.connect(Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Axis Settings"))
        self.axisModeLabel.setText(_translate("Dialog", "Axis mode"))
        self.velocityGroupBox.setTitle(_translate("Dialog", "Velocity Parameters"))
        self.referenenceWavelengthLabel.setText(_translate("Dialog", "Reference Wavelength"))
        self.redshiftGroupBox.setTitle(_translate("Dialog", "Redshift Parameters"))
        self.referenenceWavelengthLabel_2.setText(_translate("Dialog", "Amount"))

