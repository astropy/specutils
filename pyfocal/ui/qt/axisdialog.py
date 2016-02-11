# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file './source/axisdialog.ui'
#
# Created by: PyQt5 UI code generator 5.5.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(400, 300)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Dialog.sizePolicy().hasHeightForWidth())
        Dialog.setSizePolicy(sizePolicy)
        self.verticalLayout = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.formLayout_2 = QtWidgets.QFormLayout()
        self.formLayout_2.setObjectName("formLayout_2")
        self.axisModeLabel = QtWidgets.QLabel(Dialog)
        self.axisModeLabel.setObjectName("axisModeLabel")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.axisModeLabel)
        self.axisModeComboBox = QtWidgets.QComboBox(Dialog)
        self.axisModeComboBox.setObjectName("axisModeComboBox")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.axisModeComboBox)
        self.verticalLayout.addLayout(self.formLayout_2)
        self.velocityGroupBox = QtWidgets.QGroupBox(Dialog)
        self.velocityGroupBox.setObjectName("velocityGroupBox")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.velocityGroupBox)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setObjectName("formLayout")
        self.referenenceWavelengthLabel = QtWidgets.QLabel(self.velocityGroupBox)
        self.referenenceWavelengthLabel.setObjectName("referenenceWavelengthLabel")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.referenenceWavelengthLabel)
        self.referenenceWavelengthLineEdit = QtWidgets.QLineEdit(self.velocityGroupBox)
        self.referenenceWavelengthLineEdit.setObjectName("referenenceWavelengthLineEdit")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.referenenceWavelengthLineEdit)
        self.verticalLayout_2.addLayout(self.formLayout)
        self.verticalLayout.addWidget(self.velocityGroupBox)
        self.redshiftGroupBox = QtWidgets.QGroupBox(Dialog)
        self.redshiftGroupBox.setObjectName("redshiftGroupBox")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.redshiftGroupBox)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.formLayout_3 = QtWidgets.QFormLayout()
        self.formLayout_3.setObjectName("formLayout_3")
        self.referenenceWavelengthLabel_2 = QtWidgets.QLabel(self.redshiftGroupBox)
        self.referenenceWavelengthLabel_2.setObjectName("referenenceWavelengthLabel_2")
        self.formLayout_3.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.referenenceWavelengthLabel_2)
        self.referenenceWavelengthLineEdit_2 = QtWidgets.QLineEdit(self.redshiftGroupBox)
        self.referenenceWavelengthLineEdit_2.setObjectName("referenenceWavelengthLineEdit_2")
        self.formLayout_3.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.referenenceWavelengthLineEdit_2)
        self.verticalLayout_3.addLayout(self.formLayout_3)
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

