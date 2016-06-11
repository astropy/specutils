import os
from contextlib import contextmanager

QWIDGETSIZE_MAX = ((1 << 24) - 1)
ICON_PATH = os.path.abspath(os.path.join(__file__, '..', '..', 'icons'))


@contextmanager
def updates_disabled(widget):
    """Disable QWidget updates (using QWidget.setUpdatesEnabled)
    """
    old_state = widget.updatesEnabled()
    widget.setUpdatesEnabled(False)
    try:
        yield
    finally:
        widget.setUpdatesEnabled(old_state)