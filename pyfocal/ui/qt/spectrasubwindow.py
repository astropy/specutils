# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file './source/spectrasubwindow.ui'
#
# Created: Wed Jan  6 17:18:28 2016
#      by: pyside-uic 0.2.15 running on PySide 1.2.2
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_SpectraSubWindow(object):
    def setupUi(self, SpectraSubWindow):
        SpectraSubWindow.setObjectName("SpectraSubWindow")
        SpectraSubWindow.resize(800, 600)
        self.centralwidget = QtGui.QWidget(SpectraSubWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtGui.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        SpectraSubWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtGui.QStatusBar(SpectraSubWindow)
        self.statusbar.setObjectName("statusbar")
        SpectraSubWindow.setStatusBar(self.statusbar)
        self.toolBar = QtGui.QToolBar(SpectraSubWindow)
        self.toolBar.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.toolBar.setMovable(False)
        self.toolBar.setFloatable(False)
        self.toolBar.setObjectName("toolBar")
        SpectraSubWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.actionInsert_ROI = QtGui.QAction(SpectraSubWindow)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/Rectangle Stroked-50.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionInsert_ROI.setIcon(icon)
        self.actionInsert_ROI.setObjectName("actionInsert_ROI")
        self.actionInsert_Elliptical_ROI = QtGui.QAction(SpectraSubWindow)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(":/Ellipse Stroked-50.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionInsert_Elliptical_ROI.setIcon(icon1)
        self.actionInsert_Elliptical_ROI.setObjectName("actionInsert_Elliptical_ROI")
        self.actionPolygon_ROI = QtGui.QAction(SpectraSubWindow)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(":/Pentagon-50.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionPolygon_ROI.setIcon(icon2)
        self.actionPolygon_ROI.setObjectName("actionPolygon_ROI")
        self.actionGraph_Settings = QtGui.QAction(SpectraSubWindow)
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(":/Settings-50.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionGraph_Settings.setIcon(icon3)
        self.actionGraph_Settings.setObjectName("actionGraph_Settings")
        self.actionPlot_Settings = QtGui.QAction(SpectraSubWindow)
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap(":/Settings 3-50.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionPlot_Settings.setIcon(icon4)
        self.actionPlot_Settings.setObjectName("actionPlot_Settings")
        self.actionChange_Top_Axis = QtGui.QAction(SpectraSubWindow)
        self.actionChange_Top_Axis.setObjectName("actionChange_Top_Axis")
        self.actionChange_Units = QtGui.QAction(SpectraSubWindow)
        self.actionChange_Units.setObjectName("actionChange_Units")
        self.toolBar.addAction(self.actionInsert_ROI)
        self.toolBar.addAction(self.actionInsert_Elliptical_ROI)
        self.toolBar.addAction(self.actionPolygon_ROI)
        self.toolBar.addSeparator()
        self.toolBar.addAction(self.actionGraph_Settings)
        self.toolBar.addAction(self.actionPlot_Settings)

        self.retranslateUi(SpectraSubWindow)
        QtCore.QMetaObject.connectSlotsByName(SpectraSubWindow)

    def retranslateUi(self, SpectraSubWindow):
        SpectraSubWindow.setWindowTitle(QtGui.QApplication.translate("SpectraSubWindow", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.toolBar.setWindowTitle(QtGui.QApplication.translate("SpectraSubWindow", "toolBar", None, QtGui.QApplication.UnicodeUTF8))
        self.actionInsert_ROI.setText(QtGui.QApplication.translate("SpectraSubWindow", "Rectangle ROI", None, QtGui.QApplication.UnicodeUTF8))
        self.actionInsert_ROI.setToolTip(QtGui.QApplication.translate("SpectraSubWindow", "Insert a rectangular ROI selection box", None, QtGui.QApplication.UnicodeUTF8))
        self.actionInsert_Elliptical_ROI.setText(QtGui.QApplication.translate("SpectraSubWindow", "Ellipse ROI", None, QtGui.QApplication.UnicodeUTF8))
        self.actionInsert_Elliptical_ROI.setToolTip(QtGui.QApplication.translate("SpectraSubWindow", "Insert an elliptical ROI selection box", None, QtGui.QApplication.UnicodeUTF8))
        self.actionPolygon_ROI.setText(QtGui.QApplication.translate("SpectraSubWindow", "Polygon ROI", None, QtGui.QApplication.UnicodeUTF8))
        self.actionPolygon_ROI.setToolTip(QtGui.QApplication.translate("SpectraSubWindow", "Insert a polygonal ROI selection box", None, QtGui.QApplication.UnicodeUTF8))
        self.actionGraph_Settings.setText(QtGui.QApplication.translate("SpectraSubWindow", "Graph Settings", None, QtGui.QApplication.UnicodeUTF8))
        self.actionGraph_Settings.setToolTip(QtGui.QApplication.translate("SpectraSubWindow", "Edit graph settings", None, QtGui.QApplication.UnicodeUTF8))
        self.actionPlot_Settings.setText(QtGui.QApplication.translate("SpectraSubWindow", "Plot Settings", None, QtGui.QApplication.UnicodeUTF8))
        self.actionPlot_Settings.setToolTip(QtGui.QApplication.translate("SpectraSubWindow", "Edit plot settings", None, QtGui.QApplication.UnicodeUTF8))
        self.actionChange_Top_Axis.setText(QtGui.QApplication.translate("SpectraSubWindow", "Change Top Axis", None, QtGui.QApplication.UnicodeUTF8))
        self.actionChange_Units.setText(QtGui.QApplication.translate("SpectraSubWindow", "Change Units", None, QtGui.QApplication.UnicodeUTF8))

import icon_resource_rc
