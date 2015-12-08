"""
Pyfocal front-end GUI access point.

This script will start the GUI.
"""
from __future__ import absolute_import, division, print_function

import sys
from qtpy.QtGui import *

# import qdarkstyle
from pyfocal.ui.viewer import Viewer
from pyfocal.ui.controller import Controller


class App(object):
    def __init__(self):
        super(App, self).__init__()
        self.viewer = Viewer()
        self.controller = Controller(self.viewer)
        self.controller.create_sub_window()


def main():
    qapp = QApplication(sys.argv)

    # setup stylesheet
    # qapp.setStyleSheet(qdarkstyle.load_stylesheet())

    app = App()
    app.viewer.show()
    sys.exit(qapp.exec_())


if __name__ == '__main__':
    main()

    # # Start the plotting server
    # import bokeh.server
    # import bokeh.server.start
    #
    # try:
    #     bokeh.server.run()
    # except KeyboardInterrupt:
    #     bokeh.server.start.stop()
    #     print("Shutting down bokeh-server ...")

