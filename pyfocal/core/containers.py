from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import random

from ..third_party.qtpy.QtGui import *

from astropy.units import Unit, Quantity
import pyqtgraph as pg


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

        r, g, b = (random.randint(10, 25) * 10, random.randint(10, 25) * 10,
                   random.randint(10, 25) * 10)

        rand_pen = pg.mkPen(QColor(r, g, b, 255))

        _pen = pg.mkPen(pen) if pen is not None else rand_pen

        _inactive_pen = pg.mkPen(QColor(_pen.color().red(),
                                        _pen.color().green(),
                                        _pen.color().blue(),
                                        50))

        _err_pen = err_pen if err_pen is not None else pg.mkPen(
            color=(100, 100, 100, 50))

        self._pen_stash = {'pen_on': pg.mkPen(_pen),
                           'pen_inactive': pg.mkPen(_inactive_pen),
                           'pen_off': pg.mkPen(None),
                           'error_pen_on': _err_pen,
                           'error_pen_off': pg.mkPen(None)}
        self._visibility_state = [True, True, False]
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
            pen_show, _, _ = self._visibility_state

        error_pen_show = error_pen_show if pen_show else False

        if pen_show:
            self.plot.setPen(self._pen_stash['pen_on'])

            if inactive:
                self.plot.setPen(self._pen_stash['pen_inactive'])
        else:
            self.plot.setPen(self._pen_stash['pen_off'])

        if error_pen_show:
            if self.error is not None:
                self.error.setOpts(pen=self._pen_stash['error_pen_on'])
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
        _pen = self._pen_stash['pen_on']
        _inactive_pen = pg.mkPen(QColor(_pen.color().red(),
                                        _pen.color().green(),
                                        _pen.color().blue(),
                                        50))
        self._pen_stash['pen_inactive'] = _inactive_pen
        self.set_visibility(*self._visibility_state)

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
