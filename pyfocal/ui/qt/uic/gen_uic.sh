#!/bin/bash

pyuic4 ./source/mainwindow.ui -o ../mainwindow.py
pyuic4 ./source/spectrasubwindow.ui -o ../spectrasubwindow.py
pyrcc4 ./source/icon_resource.qrc -o ../icon_resource_rc.py