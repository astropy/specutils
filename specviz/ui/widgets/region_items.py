from qtpy import QtWidgets
from qtpy import QtCore

import pyqtgraph as pg


class LinearRegionItem(pg.LinearRegionItem):
    """
    Linear Region Item
    """
    sigHoverEvent = QtCore.Signal(object)
    sigRemoveRequested = QtCore.Signal(object)
    sigClicked = QtCore.Signal(object, object)

    def __init__(self, removable=True, *args, **kwargs):
        super(LinearRegionItem, self).__init__(*args, **kwargs)

        self.menu = None
        self.removable = removable
        self.deletable = removable

    def hoverEvent(self, ev):
        if self.movable and (not ev.isExit()) and ev.acceptDrags(QtCore.Qt.LeftButton):
            self.sigHoverEvent.emit(self)
            self.setMouseHover(True)
        else:
            self.setMouseHover(False)

    def contextMenuEnabled(self):
        return self.removable

    def raiseContextMenu(self, ev):
        if not self.contextMenuEnabled():
            return
        menu = self.getMenu()
        menu = self.scene().addParentContextMenus(self, menu, ev)
        pos = ev.screenPos()
        menu.popup(QtCore.QPoint(pos.x(), pos.y()))

    def getMenu(self):
        if self.menu is None:
            self.menu = QtWidgets.QMenu()
            self.menu.setTitle("ROI")
            remAct = QtWidgets.QAction("Remove ROI", self.menu)
            remAct.triggered.connect(self.removeClicked)
            self.menu.addAction(remAct)
            self.menu.remAct = remAct
        return self.menu

    def removeClicked(self):
        ## Send remove event only after we have exited the menu event handler
        QtCore.QTimer.singleShot(0, lambda: self.sigRemoveRequested.emit(self))

    def mouseClickEvent(self, ev):
        if self.moving and ev.button() == QtCore.Qt.RightButton:
            ev.accept()
            for i, l in enumerate(self.lines):
                l.setPos(self.startPositions[i])
            self.moving = False
            self.sigRegionChanged.emit(self)
            self.sigRegionChangeFinished.emit(self)
        elif int(ev.button() & self.acceptedMouseButtons()) > 0:
            ev.accept()
            if ev.button() == QtCore.Qt.RightButton and self.deletable:
                self.raiseContextMenu(ev)
            self.sigClicked.emit(self, ev)
        else:
            ev.ignore()
