"""SpecViz front-end GUI access point.
This script will start the GUI.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# STDLIB
import signal
import sys
import warnings
import os

# THIRD-PARTY
from astropy.utils.exceptions import AstropyUserWarning

# LOCAL
from .third_party.qtpy.QtWidgets import QApplication
from .third_party.qtpy.QtGui import QIcon
from .third_party.qtpy.QtCore import QTimer
from .ui.viewer import Viewer
from .ui.widgets.utils import ICON_PATH


class App(object):
    def __init__(self, argv):
        super(App, self).__init__()
        self.viewer = Viewer()

        if len(argv) > 1:
            self.controller.read_file(sys.argv[1])

def setup():
    qapp = QApplication(sys.argv)
    # qapp.setGraphicsSystem('native')
    qapp.setWindowIcon(QIcon(os.path.join(ICON_PATH, 'application',
                                          'icon.png')))

    #http://stackoverflow.com/questions/4938723/what-is-the-correct-way-to-make-my-pyqt-application-quit-when-killed-from-the-co
    timer = QTimer()
    timer.start(500)  # You may change this if you wish.
    timer.timeout.connect(lambda: None)  # Let the interpreter run each 500 ms.

    app = App(sys.argv)
    app.viewer.main_window.show()

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
    signal.signal(signal.SIGINT, sigint_handler)
    qapp, app = setup()
    sys.exit(qapp.exec_())


def sigint_handler(*args):
    """Handler for the SIGINT signal."""
    warnings.warn('KeyboardInterrupt caught; specviz will terminate',
                  AstropyUserWarning)
    QApplication.quit()


if __name__ == '__main__':
    main()
