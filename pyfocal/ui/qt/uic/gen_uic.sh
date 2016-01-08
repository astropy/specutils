#!/bin/bash

pyside-uic ./source/mainwindow.ui -o ../mainwindow.py
pyside-uic ./source/spectrasubwindow.ui -o ../spectrasubwindow.py
pyside-rcc ./source/icon_resource.qrc -o ../icon_resource_rc.py