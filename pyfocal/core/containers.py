from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .comms import EventNode
from astropy.units import Unit, Quantity
import pyqtgraph as pg
import random


class PlotContainer(object):
    def __init__(self, layer, plot=None, visible=True, style='line',
                 pen=None, err_pen=None):
        self._layer = layer
        self.style = style
        self._plot = plot
        self.error = None
        self._plot_units = None

        if self._plot is not None:
            self.change_units(self._layer.units[0], self._layer.units[1])

        rand_pen_color = pg.mkPen(
            color=(random.randint(0, 25) * 10,
                   random.randint(0, 25) * 10,
                   random.randint(0, 25) * 10,
                   255))

        print(rand_pen_color)

        _pen = pen if pen is not None else rand_pen_color
        _err_pen = err_pen if err_pen is not None else pg.mkPen(color=(0, 0, 0, 50))
        self._pen_stash = {'pen_on': _pen,
                           'pen_inactive': pg.mkPen(color=(0, 0, 0, 127)),
                           'pen_off': pg.mkPen(None),
                           'error_pen_on': _err_pen,
                           'error_pen_off': pg.mkPen(None)}
        self._visibility_state = [True, False, True]
        self.set_visibility(*self._visibility_state, override=True)

    def change_units(self, x, y=None, z=None):
        self._plot_units = (
            x or self._plot_units[0] or self.layer.layer_units[0],
            y or self._plot_units[1] or self.layer.layer_units[1],
            z)

        self.update()

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
    def data(self):
        return self.layer.data.to(self._plot_units[1])

    @property
    def dispersion(self):
        return self.layer.dispersion.to(self._plot_units[0])

    @property
    def uncertainty(self):
        return Quantity(self.layer.uncertainty.array,
                        unit=self.layer.units[1]).to(self._plot_units[1])

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
        self._pen_stash['pen_on'] = pg.mkPen(pen)
        self.plot.setPen(pg.mkPen(pen))

    @property
    def error_pen(self):
        return self._pen_stash['error_pen_on']

    @error_pen.setter
    def error_pen(self, pen):
        self._pen_stash['error_pen_on'] = pg.mkPen(pen)

        if self.error is not None:
            self.error.setOpts(pen=pg.mkPen(pen))

    def update(self, autoscale=False):
        self._plot.setData(self.dispersion.value,
                           self.data.value)

        if self.error is not None:
            self.error.setData(
                    x=self.dispersion.value,
                    y=self.data.value,
                    height=self.uncertainty.value)
