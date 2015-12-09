import pyqtgraph as pg


class Profile(pg.PlotWidget):
    def __init__(self, data=None, *args, **kwargs):
        super(Profile, self).__init__(*args, **kwargs)
        self._plot_item = self.getPlotItem()

        if data is not None:
            self._plot_item.plot(data.data)

