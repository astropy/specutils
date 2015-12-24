from __future__ import absolute_import, division, print_function

from .base_plot import BasePlot
from pyfocal.core.containers import PlotContainer


class Plot(BasePlot):
    """
    One-dimensional representation of a set of data.
    """
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
        plot_container = PlotContainer(layer=layer, unit=unit or layer.unit,
                                       visible=visible, style=style, pen=pen)
        plot = self._plot_item.plot(plot_container.data)
        plot_container.plot = plot

        self._containers.append(plot_container)

    def change_unit(self, new_unit):
        for plot_container in self._containers:
            plot_container.change_unit(new_unit)

        self._plot_item.setLabels(bottom=str(self._containers[0].unit))
