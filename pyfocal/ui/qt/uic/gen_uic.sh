#!/bin/bash

case $1 in
    pyside)
    uic=pyside-uic
    rcc=pyside-rcc
    ;;
    pyqt4)
    uic=pyuic4
    rcc=pyrcc4
    ;;
    pyqt5)
    uic=pyuic5
    rcc=pyrcc5
    ;;
esac

$uic ./source/mainwindow.ui -o ../mainwindow.py
$uic ./source/plotsubwindow.ui -o ../plotsubwindow.py
$rcc ./source/icon_resource.qrc -o ../icon_resource_rc.py
