#-------------------------------------------------
#
# Project created by QtCreator 2015-09-30T22:55:57
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = specviz
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp

HEADERS  += mainwindow.h

FORMS    += mainwindow.ui \
    plotsubwindow.ui \
    axisdialog.ui \
    unit_change_dialog.ui \
    layer_arithmetic_dialog.ui

RESOURCES += \
    icon_resource.qrc

DISTFILES += \
    plot_sub_window_plugin.py
