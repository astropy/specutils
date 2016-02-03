import pyqtgraph as pg
import numpy as np
from functools import reduce

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
pg.setConfigOptions(antialias=False)


class Plot(pg.PlotWidget):
    def __init__(self, data=None, parent=None, *args, **kwargs):
        super(Plot, self).__init__(*args, **kwargs)
        self._parent = parent

        # Cache plot item
        self._plot_item = self.getPlotItem()

        # Create roi container
        self._rois = []

        self._setup_connections()

        # Add data
        if data is not None:
            self.add_data(data)

    def _setup_connections(self):
        if self._parent is not None:
            act_insert_roi = self._parent.action("actionInsert_ROI")
            act_insert_roi.triggered.connect(self.add_roi)

    def add_roi(self):
        view_range = self.viewRange()
        x_len = (view_range[0][1] - view_range[0][0]) * 0.5
        y_len = (view_range[1][1] - view_range[1][0]) * 0.5
        x_pos = x_len + view_range[0][0]
        y_pos = y_len + view_range[1][0]

        def remove():
            self.removeItem(roi)
            self._rois.remove(roi)

        roi = pg.RectROI([x_pos, y_pos], [x_len * 0.5, y_len * 0.5],
                         sideScalers=True, removable=True, pen='k',
                         hoverPen='r', handlePen='k')
        self._rois.append(roi)
        self.addItem(roi)

        # Connect the remove functionality
        roi.sigRemoveRequested.connect(remove)

        # Let everyone know, ROI is ready for use.
        roi.sigRegionChangeFinished.emit(roi)

    def add_data(self, data):
        raise NotImplemented()

    def update(self):
        raise NotImplemented()

    def get_roi_mask(self, container):
        mask_holder = []

        for roi in self._rois:
            roi_shape = roi.parentBounds()
            x1, y1, x2, y2 = roi_shape.getCoords()

            mask_holder.append((container.layer.dispersion.value >= x1) &
                               (container.layer.dispersion.value <= x2) &
                               (container.layer.data.value >= y1) &
                               (container.layer.data.value <= y2))

        if len(mask_holder) == 0:
            mask_holder.append(np.ones(shape=container.dispersion.value.shape,
                                       dtype=bool))

        # mask = np.logical_not(reduce(np.logical_or, mask_holder))
        mask = reduce(np.logical_or, mask_holder)
        return mask
