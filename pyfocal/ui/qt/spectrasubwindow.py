# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file './source/spectrasubwindow.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_SpectraSubWindow(object):
    def setupUi(self, SpectraSubWindow):
        SpectraSubWindow.setObjectName(_fromUtf8("SpectraSubWindow"))
        SpectraSubWindow.resize(800, 600)
        self.centralwidget = QtGui.QWidget(SpectraSubWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.gridLayout = QtGui.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        SpectraSubWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtGui.QStatusBar(SpectraSubWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        SpectraSubWindow.setStatusBar(self.statusbar)
        self.toolBar = QtGui.QToolBar(SpectraSubWindow)
        self.toolBar.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.toolBar.setMovable(False)
        self.toolBar.setFloatable(False)
        self.toolBar.setObjectName(_fromUtf8("toolBar"))
        SpectraSubWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.actionInsert_ROI = QtGui.QAction(SpectraSubWindow)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/Rectangle Stroked-50.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionInsert_ROI.setIcon(icon)
        self.actionInsert_ROI.setObjectName(_fromUtf8("actionInsert_ROI"))
        self.actionInsert_Elliptical_ROI = QtGui.QAction(SpectraSubWindow)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(_fromUtf8(":/Ellipse Stroked-50.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionInsert_Elliptical_ROI.setIcon(icon1)
        self.actionInsert_Elliptical_ROI.setObjectName(_fromUtf8("actionInsert_Elliptical_ROI"))
        self.actionPolygon_ROI = QtGui.QAction(SpectraSubWindow)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(_fromUtf8(":/Pentagon-50.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionPolygon_ROI.setIcon(icon2)
        self.actionPolygon_ROI.setObjectName(_fromUtf8("actionPolygon_ROI"))
        self.actionGraph_Settings = QtGui.QAction(SpectraSubWindow)
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(_fromUtf8(":/Settings-50.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionGraph_Settings.setIcon(icon3)
        self.actionGraph_Settings.setObjectName(_fromUtf8("actionGraph_Settings"))
        self.actionPlot_Settings = QtGui.QAction(SpectraSubWindow)
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap(_fromUtf8(":/Settings 3-50.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionPlot_Settings.setIcon(icon4)
        self.actionPlot_Settings.setObjectName(_fromUtf8("actionPlot_Settings"))
        self.actionChange_Top_Axis = QtGui.QAction(SpectraSubWindow)
        self.actionChange_Top_Axis.setObjectName(_fromUtf8("actionChange_Top_Axis"))
        self.actionChange_Units = QtGui.QAction(SpectraSubWindow)
        self.actionChange_Units.setObjectName(_fromUtf8("actionChange_Units"))
        self.toolBar.addAction(self.actionInsert_ROI)
        self.toolBar.addAction(self.actionInsert_Elliptical_ROI)
        self.toolBar.addAction(self.actionPolygon_ROI)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionGraph_Settings)
        self.toolBar.addAction(self.actionPlot_Settings)

        self.retranslateUi(SpectraSubWindow)
        QtCore.QMetaObject.connectSlotsByName(SpectraSubWindow)

    def retranslateUi(self, SpectraSubWindow):
        SpectraSubWindow.setWindowTitle(_translate("SpectraSubWindow", "MainWindow", None))
        self.toolBar.setWindowTitle(_translate("SpectraSubWindow", "toolBar", None))
        self.actionInsert_ROI.setText(_translate("SpectraSubWindow", "Rectangle ROI", None))
        self.actionInsert_ROI.setToolTip(_translate("SpectraSubWindow", "Insert a rectangular ROI selection box", None))
        self.actionInsert_Elliptical_ROI.setText(_translate("SpectraSubWindow", "Ellipse ROI", None))
        self.actionInsert_Elliptical_ROI.setToolTip(_translate("SpectraSubWindow", "Insert an elliptical ROI selection box", None))
        self.actionPolygon_ROI.setText(_translate("SpectraSubWindow", "Polygon ROI", None))
        self.actionPolygon_ROI.setToolTip(_translate("SpectraSubWindow", "Insert a polygonal ROI selection box", None))
        self.actionGraph_Settings.setText(_translate("SpectraSubWindow", "Graph Settings", None))
        self.actionGraph_Settings.setToolTip(_translate("SpectraSubWindow", "Edit graph settings", None))
        self.actionPlot_Settings.setText(_translate("SpectraSubWindow", "Plot Settings", None))
        self.actionPlot_Settings.setToolTip(_translate("SpectraSubWindow", "Edit plot settings", None))
        self.actionChange_Top_Axis.setText(_translate("SpectraSubWindow", "Change Top Axis", None))
        self.actionChange_Units.setText(_translate("SpectraSubWindow", "Change Units", None))

import icon_resource_rc
