# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file './source/plotsubwindow.ui'
#
# Created by: ...third_party.qtpy UI code generator 5.5.1
#
# WARNING! All changes made in this file will be lost!
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...third_party.qtpy import QtCore, QtGui, QtWidgets

class Ui_SpectraSubWindow(object):
    def setupUi(self, SpectraSubWindow):
        SpectraSubWindow.setObjectName("SpectraSubWindow")
        SpectraSubWindow.resize(800, 600)
        SpectraSubWindow.setWindowTitle("")
        self.centralwidget = QtWidgets.QWidget(SpectraSubWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        SpectraSubWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(SpectraSubWindow)
        self.statusbar.setObjectName("statusbar")
        SpectraSubWindow.setStatusBar(self.statusbar)
        self.toolBar = QtWidgets.QToolBar(SpectraSubWindow)
        self.toolBar.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.toolBar.setMovable(False)
        self.toolBar.setFloatable(False)
        self.toolBar.setObjectName("toolBar")
        SpectraSubWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.actionInsert_ROI = QtWidgets.QAction(SpectraSubWindow)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/Rectangle Stroked-50.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionInsert_ROI.setIcon(icon)
        self.actionInsert_ROI.setObjectName("actionInsert_ROI")
        self.actionInsert_Elliptical_ROI = QtWidgets.QAction(SpectraSubWindow)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(":/Ellipse Stroked-50.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionInsert_Elliptical_ROI.setIcon(icon1)
        self.actionInsert_Elliptical_ROI.setObjectName("actionInsert_Elliptical_ROI")
        self.actionPolygon_ROI = QtWidgets.QAction(SpectraSubWindow)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(":/Pentagon-50.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionPolygon_ROI.setIcon(icon2)
        self.actionPolygon_ROI.setObjectName("actionPolygon_ROI")
        self.actionGraph_Settings = QtWidgets.QAction(SpectraSubWindow)
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(":/Settings-50.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionGraph_Settings.setIcon(icon3)
        self.actionGraph_Settings.setObjectName("actionGraph_Settings")
        self.actionPlot_Settings = QtWidgets.QAction(SpectraSubWindow)
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap(":/Settings 3-50.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionPlot_Settings.setIcon(icon4)
        self.actionPlot_Settings.setObjectName("actionPlot_Settings")
        self.actionChange_Top_Axis = QtWidgets.QAction(SpectraSubWindow)
        self.actionChange_Top_Axis.setObjectName("actionChange_Top_Axis")
        self.actionChange_Units = QtWidgets.QAction(SpectraSubWindow)
        self.actionChange_Units.setObjectName("actionChange_Units")
        self.toolBar.addAction(self.actionInsert_ROI)
        self.toolBar.addSeparator()

        self.retranslateUi(SpectraSubWindow)
        QtCore.QMetaObject.connectSlotsByName(SpectraSubWindow)

    def retranslateUi(self, SpectraSubWindow):
        _translate = QtCore.QCoreApplication.translate
        self.toolBar.setWindowTitle(_translate("SpectraSubWindow", "toolBar"))
        self.actionInsert_ROI.setText(_translate("SpectraSubWindow", "Rectangle ROI"))
        self.actionInsert_ROI.setToolTip(_translate("SpectraSubWindow", "Insert a rectangular ROI selection box"))
        self.actionInsert_Elliptical_ROI.setText(_translate("SpectraSubWindow", "Ellipse ROI"))
        self.actionInsert_Elliptical_ROI.setToolTip(_translate("SpectraSubWindow", "Insert an elliptical ROI selection box"))
        self.actionPolygon_ROI.setText(_translate("SpectraSubWindow", "Polygon ROI"))
        self.actionPolygon_ROI.setToolTip(_translate("SpectraSubWindow", "Insert a polygonal ROI selection box"))
        self.actionGraph_Settings.setText(_translate("SpectraSubWindow", "Graph Settings"))
        self.actionGraph_Settings.setToolTip(_translate("SpectraSubWindow", "Edit graph settings"))
        self.actionPlot_Settings.setText(_translate("SpectraSubWindow", "Plot Settings"))
        self.actionPlot_Settings.setToolTip(_translate("SpectraSubWindow", "Edit plot settings"))
        self.actionChange_Top_Axis.setText(_translate("SpectraSubWindow", "Change Top Axis"))
        self.actionChange_Units.setText(_translate("SpectraSubWindow", "Change Units"))

from . import icon_resource_rc
