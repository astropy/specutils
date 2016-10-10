from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..third_party.qtpy.QtGui import *

from astropy.units import spectral_density, spectral
import pyqtgraph as pg
from itertools import cycle
import logging

AVAILABLE_COLORS = cycle([(0, 0, 0),
                          (0.86, 0.37119999999999997, 0.33999999999999997),
                          (0.86, 0.68320000000000003, 0.33999999999999997),
                          (0.72479999999999989, 0.86, 0.33999999999999997),
                          (0.41279999999999994, 0.86, 0.33999999999999997),
                          (0.33999999999999997, 0.86, 0.57920000000000016),
                          (0.33999999999999997, 0.82879999999999987, 0.86),
                          (0.33999999999999997, 0.51679999999999948, 0.86),
                          (0.47520000000000029, 0.33999999999999997, 0.86),
                          (0.7871999999999999, 0.33999999999999997, 0.86),
                          (0.86, 0.33999999999999997, 0.62079999999999991)])


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
        r, g, b = r * 255, g * 255, b * 255

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
            x = None

        if y is None or not self._layer.unit.is_equivalent(
                y, equivalencies=spectral_density(self.layer.dispersion)):
            logging.error("Failed to convert y-axis plot units.")
            y = None

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
        if hasattr(self.layer, '_model'):
            disp = self.layer.unmasked_dispersion.compressed().value
            data = self.layer.unmasked_data.compressed().value
            uncert = self.layer.unmasked_raw_uncertainty.compressed().value
        else:
            disp = self.layer.dispersion.compressed().value
            data = self.layer.data.compressed().value
            uncert = self.layer.raw_uncertainty.compressed().value

        self._plot.setData(disp, data)

        if self.error is not None:
            self.error.setData(x=disp, y=data, height=uncert)
