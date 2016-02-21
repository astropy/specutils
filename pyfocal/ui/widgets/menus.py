from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...third_party.qtpy.QtWidgets import *


class BaseContextMenu(QMenu):
    def __init__(self, *args, **kwargs):
        super(BaseContextMenu, self).__init__(*args, **kwargs)

        self.act_change_color = self.addAction("Change color")
        self.act_export = self.addAction("Export")
        self.act_export.setDisabled(True)

class LayerContextMenu(BaseContextMenu):
    def __init__(self, *args, **kwargs):
        super(LayerContextMenu, self).__init__(*args, **kwargs)
