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

ui_files=('mainwindow' 'plotsubwindow' 'axisdialog' 'layer_arithmetic_dialog' 'unit_change_dialog')
rc_files=('icon_resource')

for i in "${ui_files[@]}"
do
    $uic ./source/${i}.ui -o ../${i}.py;1
    sed -i '' 's/import icon_resource_rc/from . import icon_resource_rc/g; s/PyQt5/...third_party.qtpy/g' ../${i}.py;
    sed -i '' '8s/^/from __future__ import (absolute_import, division, print_function,\
                        unicode_literals)/' ../$i.py
done

for i in "${rc_files[@]}"
do
    $rcc ./source/${i}.qrc -o ../${i}_rc.py;
    sed -i '' 's/import icon_resource_rc/from . import icon_resource_rc/g; s/PyQt5/...third_party.qtpy/g' ../${i}_rc.py;
done
