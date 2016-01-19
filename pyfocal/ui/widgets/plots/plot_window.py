from .base_plot import BasePlot
from pyfocal.core.containers import PlotContainer

import pyqtgraph as pg
import random


class PlotWindow(BasePlot):
    """
    One-dimensional representation of a set of data.
    """
    inactive_color = pg.mkPen(color=(0, 0, 0, 75))
    active_color = pg.mkPen(color=(0, 0, 0, 255))

    def __init__(self, *args, **kwargs):
        super(PlotWindow, self).__init__(*args, **kwargs)

        self.highlight = self._plot_item.plot()

    def add_data(self, layer, unit=None, visible=False, style='line',
                 pen=None):
        """
        This method creates a container with the plot settings particular
        to this data set.

        Parameters
        ----------
        layer : pyfocal.core.data.Layer
            Layer object from which the plot will retrieve data.
        """
        if pen is None:
            print("Pen is not set.")
            pen = pg.mkPen()

        plot_container = PlotContainer(layer=layer, visible=visible,
                                       style=style, pen=pen)
        plot = self._plot_item.plot(
                plot_container.dispersion.value,
                plot_container.data.value)
        plot_container.plot = plot

        self._containers.append(plot_container)
        self.set_labels()

    def change_unit(self, new_unit):
        for plot_container in self._containers:
            plot_container.change_unit(new_unit)

    def set_labels(self):
        self._plot_item.setLabels(bottom=str(self._containers[0].units[0]))

    def set_active_plot(self, layer):
        for container in self._containers:
            if container.layer == layer:
                container.set_pen(self.active_color)
            else:
                container.set_pen(self.inactive_color)
