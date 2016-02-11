from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .events import EventHook
from astropy.units import Unit, Quantity


class PlotContainer(object):
    def __init__(self, layer, plot=None, visible=True, style='line',
                 pen=None, err_pen=None):
        self._layer = layer
        self.visible = visible
        self.style = style
        self._pen = pen
        self._err_pen = err_pen
        self._plot = plot
        self.error = None

        self.on_unit_change = EventHook()
        self.on_visibility_change = EventHook()
        self.on_pen_change = EventHook()

    def change_unit(self, x, y=None, z=None):
        self.layer.units = (x, y or self.layer.layer_units[1],
                            z or self.layer.layer_units[2])

    def change_visible(self, visible_state):
        self.visible = visible_state

    @property
    def plot(self):
        return self._plot

    @plot.setter
    def plot(self, plot_item):
        self._plot = plot_item
        # self._plot.setPen(self.pen)

    @property
    def layer(self):
        return self._layer

    @property
    def pen(self):
        return self._pen

    @pen.setter
    def pen(self, pen):
        self._pen = pen
        self.plot.setPen(self._pen)

    @property
    def error_pen(self):
        return self._err_pen

    @error_pen.setter
    def error_pen(self, pen):
        self._err_pen = pen

        if self.error is not None:
            # self.error.setPen(self._err_pen)
            # self.error.setBrush(self._err_pen)
            self.error.setOpts(pen=pen)

    def update(self):
        print("PlotContainer.update: {}".format(type(self._layer)))
        self._plot.setData(self._layer.dispersion.value,
                           self._layer.data.value)

        if self.error is not None:
            self.error.setData(
                    x=self._layer.dispersion.value,
                    y=self._layer.data.value,
                    height=self._layer.uncertainty.array)
