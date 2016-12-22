"""
Base contextual menu classes
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from qtpy.QtWidgets import *


class BaseContextMenu(QMenu):
    """
    Base class for all context menus
    """
    def __init__(self, *args, **kwargs):
        super(BaseContextMenu, self).__init__(*args, **kwargs)


class LayerContextMenu(BaseContextMenu):
    """
    Base class for Layer-based contextual menus
    """
    def __init__(self, *args, **kwargs):
        super(LayerContextMenu, self).__init__(*args, **kwargs)

        self.act_change_color = self.addAction("Change color")
        self.act_export = self.addAction("Export")
        self.act_export.setDisabled(True)


class ModelContextMenu(BaseContextMenu):
    """
    Base class for all Model-based contextual menus
    """
    def __init__(self, *args, **kwargs):
        super(ModelContextMenu, self).__init__(*args, **kwargs)

        self.act_remove = self.addAction("Remove")
