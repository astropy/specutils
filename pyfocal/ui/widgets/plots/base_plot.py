import pyqtgraph as pg

import numpy as np

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
pg.setConfigOptions(antialias=False)


class BasePlot(pg.PlotWidget):
    def __init__(self, data=None, parent=None, *args, **kwargs):
        super(BasePlot, self).__init__(*args, **kwargs)
        self.parent = parent
        self._containers = []

        # Cache plot item
        self._plot_item = self.getPlotItem()

        # Create roi container
        self._rois = []

        self._setup_connections()

        # Add data
        if data is not None:
            self.add_data(data)

    def _setup_connections(self):
        self.parent.actionInsert_ROI.triggered.connect(self.add_roi)

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
                         sideScalers=True, removable=True, pen='k')
        self._rois.append(roi)
        self.addItem(roi)

        # Connect the remove functionality
        roi.sigRemoveRequested.connect(remove)

        # Let everyone know, ROI is ready for use.
        roi.sigRegionChangeFinished.emit(self)

    def add_data(self, data):
        raise NotImplemented()

    def update(self):
        raise NotImplemented()

    def get_roi_mask(self, layer):
        data = self.get_container(layer).data

        mask_holder = []

        for roi in self._rois:
            roi_shape = roi.parentBounds()
            x1, y1, x2, y2 = roi_shape.getCoords()

            mask_holder.append((x_data.value >= x1) & (x_data.value <= x2) &
                               (y_data.value >= y1) & (y_data.value <= y2))
        else:
            mask_holder.append(np.ones(shape=x_data.value.shape, dtype=bool))

        # mask = np.logical_not(reduce(np.logical_or, mask_holder))
        mask = reduce(np.logical_or, mask_holder)

        return mask

    def get_container(self, layer):
        """
        Retrieve the container of the specified layer.

        Parameters
        ----------
        layer : pyfocal.core.data.Layer
            The layer for which to get the container.
        """
        for container in self._containers:
            if container.layer == layer:
                return layer

        return None