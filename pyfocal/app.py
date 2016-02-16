"""
Pyfocal front-end GUI access point.

This script will start the GUI.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys

from .third_party.qtpy.QtWidgets import *

from .ui.viewer import Viewer
from .ui.controller import Controller


class App(object):
    def __init__(self, argv):
        super(App, self).__init__()
        self.viewer = Viewer()
        self.controller = Controller(self.viewer)

        if len(argv) > 1:
            self.controller.read_file(sys.argv[1])

def setup():
    qapp = QApplication(sys.argv)
    # qapp.setGraphicsSystem('native')

    app = App(sys.argv)
    app.viewer.show()

    return qapp, app


def embed():
    """
    Used when launching the application within a shell, and the application
    namespace is still needed.
    """
    qapp, app = setup()
    qapp.exec_()

    return app


def main():
    """
    Used when launching the application as standalone.
    """
    qapp, app = setup()
    sys.exit(qapp.exec_())


if __name__ == '__main__':
    main()
