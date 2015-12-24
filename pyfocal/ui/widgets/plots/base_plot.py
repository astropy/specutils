import pyqtgraph as pg

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