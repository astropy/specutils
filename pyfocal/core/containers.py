from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .events import EventHook
from astropy.units import Unit, Quantity
import pyqtgraph as pg


class PlotContainer(object):
    def __init__(self, layer, plot=None, visible=True, style='line',
                 pen=None, err_pen=None):
        self._layer = layer
        self.style = style
        self._plot = plot
        self.error = None

        _pen = pen if pen is not None else pg.mkPen(color=(0, 0, 0, 255))
        _err_pen = err_pen if err_pen is not None else pg.mkPen(color=(0, 0, 0, 50))
        self._pen_stash = {'pen_on': _pen,
                           'pen_inactive': pg.mkPen(color=(0, 0, 0, 127)),
                           'pen_off': pg.mkPen(None),
                           'error_pen_on': _err_pen,
                           'error_pen_off': pg.mkPen(None)}
        self._visibility_state = [True, True, True]

        self.on_unit_change = EventHook()
        self.on_visibility_change = EventHook()
        self.on_pen_change = EventHook()

    def change_unit(self, x, y=None, z=None):
        self.layer.units = (x, y or self.layer.layer_units[1],
                            z or self.layer.layer_units[2])

    def set_visibility(self, pen_show, error_pen_show, inactive=True,
                       override=False):
        if override:
            self._visibility_state = [pen_show, error_pen_show, inactive]
        else:
            pen_show, _, inactive = self._visibility_state

        error_pen_show = error_pen_show if pen_show else False

        if pen_show:
            self.pen = self._pen_stash['pen_on']
        else:
            if inactive:
                self.pen = self._pen_stash['pen_inactive']
            else:
                self.plot.setPen(self._pen_stash['pen_off'])

        if error_pen_show:
            self.error_pen = self._pen_stash['error_pen_on']
        else:
            if self.error is not None:
                self.error.setOpts(pen=self._pen_stash['error_pen_off'])


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
        return self._pen_stash['pen_on']

    @pen.setter
    def pen(self, pen):
        self._pen_stash['pen_on'] = pen
        self.plot.setPen(pen)

    @property
    def error_pen(self):
        return self._pen_stash['error_pen_on']

    @error_pen.setter
    def error_pen(self, pen):
        self._pen_stash['error_pen_on'] = pen

        if self.error is not None:
            self.error.setOpts(pen=pen)

    def update(self):
        self._plot.setData(self._layer.dispersion.value,
                           self._layer.data.value)

        if self.error is not None:
            self.error.setData(
                    x=self._layer.dispersion.value,
                    y=self._layer.data.value,
                    height=self._layer.uncertainty.array)
