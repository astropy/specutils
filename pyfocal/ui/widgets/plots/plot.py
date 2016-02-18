from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ....core.events import EventHook
from ..region_items import LinearRegionItem

import pyqtgraph as pg
import numpy as np
from functools import reduce

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
pg.setConfigOptions(antialias=False)


class Plot(pg.PlotWidget):
    def __init__(self, data=None, parent=None, *args, **kwargs):
        super(Plot, self).__init__(*args, **kwargs)
        # Events
        self.on_roi_update = EventHook()
        self.on_roi_add = EventHook()
        self.on_roi_remove = EventHook()
        self.on_roi_change_start = EventHook()
        self.on_roi_change_end = EventHook()

        self._parent = parent

        # Cache plot item
        self._plot_item = self.getPlotItem()

        # Create roi container
        self._rois = []
        self._active_roi = None

        # Add data
        if data is not None:
            self.add_data(data)

        # Send a statistics update event
        self.on_roi_update.emit(None)

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

        # Let everyone know, ROI is ready for use.
        roi.sigRegionChangeFinished.emit(roi)

        # Set the active ROI as the last one interacted with
        roi.sigRegionChangeFinished.connect(self._set_active_roi)

        # Connect events
        self.on_roi_update.emit(roi)
        roi.sigRemoveRequested.connect(lambda: self.on_roi_update.emit(roi))
        roi.sigRegionChangeFinished.connect(
            lambda: self.on_roi_update.emit(roi))
        # roi.sigRegionChangeStarted.connect(
        #     lambda: self.on_roi_update.emit(roi))

    def _set_active_roi(self, roi):
        self._active_roi = roi

    def add_data(self, data):
        raise NotImplemented()

    def update(self):
        raise NotImplemented()

    def get_roi_mask(self, container):
        if container is None:
            return

        mask_holder = []

        for roi in self._rois:
            # roi_shape = roi.parentBounds()
            # x1, y1, x2, y2 = roi_shape.getCoords()
            x1, x2 = roi.getRegion()

            mask_holder.append((container.layer.dispersion.value >= x1) &
                               (container.layer.dispersion.value <= x2))

        if len(mask_holder) == 0:
            mask_holder.append(np.ones(
                shape=container.layer.dispersion.value.shape,
                dtype=bool))

        # mask = np.logical_not(reduce(np.logical_or, mask_holder))
        mask = reduce(np.logical_or, mask_holder)
        return mask
