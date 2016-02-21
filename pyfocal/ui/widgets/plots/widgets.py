from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ....core.comms import Dispatch, DispatchHandle
from ..region_items import LinearRegionItem

import pyqtgraph as pg
import numpy as np
from functools import reduce

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
pg.setConfigOptions(antialias=False)


class PlotWidget(pg.PlotWidget):
    def __init__(self, data=None, window=None, *args, **kwargs):
        super(PlotWidget, self).__init__(*args, **kwargs)
        self._window = window

        # Cache plot item
        self._plot_item = self.getPlotItem()

        # Create roi container
        self._rois = []
        self._active_roi = None

        # Add data
        if data is not None:
            self.add_data(data)

        DispatchHandle.setup(self)

    def add_roi(self):
        view_range = self.viewRange()
        x_len = (view_range[0][1] - view_range[0][0]) * 0.5
        y_len = (view_range[1][1] - view_range[1][0]) * 0.9
        x_pos = x_len * 0.5 + view_range[0][0]
        y_pos = y_len * 0.05 + view_range[1][0]

        def remove():
            self.removeItem(roi)
            self._rois.remove(roi)

        roi = LinearRegionItem(values=[x_pos, x_pos + x_len])
        self._rois.append(roi)
        self.addItem(roi)

        # Connect the remove functionality
        roi.sigRemoveRequested.connect(remove)

        # Set the active ROI as the last one interacted with
        roi.sigRegionChangeFinished.connect(self._set_active_roi)

        # Connect events
        Dispatch.on_update_roi.emit(roi=roi)
        roi.sigRemoveRequested.connect(
            lambda: Dispatch.on_update_roi.emit(roi=roi))
        roi.sigRegionChangeFinished.connect(
            lambda: Dispatch.on_update_roi.emit(roi=roi))

    def _set_active_roi(self, roi):
        self._active_roi = roi

    def add_data(self, data):
        raise NotImplemented()

    def update(self):
        raise NotImplemented()
