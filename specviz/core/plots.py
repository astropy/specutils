from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import random

from ..third_party.qtpy.QtGui import *

from astropy.units import (Unit, Quantity, spectral_density, spectral,
                           UnitConversionError)
import pyqtgraph as pg
from itertools import cycle
import logging

AVAILABLE_COLORS = cycle([(0, 0, 0), (0, 73, 73), (0, 146, 146),
                          (255, 109, 182), (255, 182, 219), (73, 0, 146),
                          (0, 109, 219), (182, 109, 255), (109, 182, 255),
                          (182, 219, 255), (146, 0, 0), (146, 73, 0),
                          (219, 209, 0), (36, 255, 36), (255, 255, 109)])


class LinePlot(object):
    def __init__(self, layer, plot=None, visible=True, style='line',
                 pen=None, err_pen=None):
        self._layer = layer
        self.style = style
        self._plot = plot
        self.error = None
        self._plot_units = (self._layer.dispersion_unit,
                            self._layer.unit,
                            None)
        self.line_width = 1

        if self._plot is not None:
            self.change_units(self._layer.dispersion_unit,
                              self._layer.unit)

        r, g, b = next(AVAILABLE_COLORS)

        rand_pen = pg.mkPen(QColor(r, g, b, 255), width=self.line_width)

        _pen = pg.mkPen(pen, width=self.line_width) if pen is not None else rand_pen

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

    @staticmethod
    def from_layer(layer, **kwargs):
        plot_data_item = pg.PlotDataItem(layer.dispersion, layer.data)

        plot_container = LinePlot(layer=layer, plot=plot_data_item, **kwargs)

        if plot_container.layer.raw_uncertainty is not None:
            # err_top = pg.PlotDataItem(
            #     plot_container.layer.dispersion.value,
            #     plot_container.layer.data.value +
            #     plot_container.layer.uncertainty.array * 0.5)
            # err_btm = pg.PlotDataItem(
            #     plot_container.layer.dispersion.value,
            #     plot_container.layer.data.value -
            #     plot_container.layer.uncertainty.array * 0.5)
            #
            # plot_error_item = pg.FillBetweenItem(err_top, err_btm, 'r')
            plot_error_item = pg.ErrorBarItem(
                x=plot_container.layer.dispersion.compressed().value,
                y=plot_container.layer.data.compressed().value,
                height=plot_container.layer.raw_uncertainty.compressed().value,
            )
            plot_container.error = plot_error_item

        return plot_container

    def change_units(self, x, y=None, z=None):
        if x is None or not self._layer.dispersion_unit.is_equivalent(
                x, equivalencies=spectral()):
            logging.error("Failed to convert x-axis plot units. {} to {"
                          "}".format(self._layer.dispersion_unit, x))
            x = self._plot_units[0] or self._layer.units[0]

        if y is None or not self._layer.unit.is_equivalent(
                y, equivalencies=spectral_density(self.dispersion)):
            logging.error("Failed to convert y-axis plot units.")
            y = self._plot_units[1] or self._layer.units[1]

        self._layer.set_units(x, y)
        self._plot_units = (x, y, z)
        self.update()

    def set_visibility(self, pen_show, error_pen_show, inactive=True,
                       override=False):
        if override:
            self._visibility_state = [pen_show, error_pen_show, inactive]
        else:
            pen_show, _, _ = self._visibility_state

        error_pen_show = error_pen_show if pen_show else False

        if pen_show:
            self._plot.setPen(self._pen_stash['pen_on'])

            if inactive:
                self._plot.setPen(self._pen_stash['pen_inactive'])
        else:
            self._plot.setPen(self._pen_stash['pen_off'])

        if error_pen_show:
            if self.error is not None:
                self.error.setOpts(pen=self._pen_stash['error_pen_on'])
        else:
            if self.error is not None:
                self.error.setOpts(pen=self._pen_stash['error_pen_off'])

    @property
    def data(self):
        return self.layer.data

    @property
    def dispersion(self):
        return self.layer.dispersion

    @property
    def uncertainty(self):
        return self.layer.raw_uncertainty

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
        self._pen_stash['pen_on'] = pg.mkPen(pen, width=self.line_width)
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
        self._plot.setData(self.dispersion.compressed().value,
                           self.data.compressed().value)

        if self.error is not None:
            self.error.setData(
                    x=self.dispersion.compressed().value,
                    y=self.data.compressed().value,
                    height=self.uncertainty.compressed().value)
