import pyqtgraph as pg


class Profile(pg.PlotWidget):
    def __init__(self, data=None, **kwargs):
        self._plot_item = self.getPlotItem()

        if data is not None:
            self._plot_item.plot(data)

