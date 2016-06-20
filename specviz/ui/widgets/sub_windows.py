from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
import logging
from functools import reduce

import numpy as np
import pyqtgraph as pg

from astropy.units import Quantity

from ...third_party.qtpy.QtWidgets import *
from ...third_party.qtpy.QtCore import *

from ...core.comms import Dispatch, DispatchHandle
from ...core.linelist import LineList
from ...core.plots import LinePlot
from ...core.annotation import LineIDMarker
from .axes import DynamicAxisItem
from .region_items import LinearRegionItem

from .dialogs import LineListsWindow


pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
pg.setConfigOptions(antialias=False)


class UiPlotSubWindow(QMainWindow):
    def __init__(self, *args, **kwargs):
        super(UiPlotSubWindow, self).__init__(*args, **kwargs)

        self.vertical_layout = QVBoxLayout()
        self.horizontal_layout = QHBoxLayout()
        self.vertical_layout.setContentsMargins(0, 0, 0, 0)
        self.vertical_layout.setSpacing(2)

        # X range
        self.label_x_range = QLabel()
        self.label_x_range.setText("X Range")
        self.line_edit_min_x_range = QLineEdit()
        self.line_edit_max_x_range = QLineEdit()

        self.layout_x_range = QHBoxLayout()
        self.layout_x_range.addWidget(self.label_x_range)
        self.layout_x_range.addWidget(self.line_edit_min_x_range)
        self.layout_x_range.addWidget(self.line_edit_max_x_range)

        # Y range
        self.label_y_range = QLabel()
        self.label_y_range.setText("Y Range")
        self.line_edit_min_y_range = QLineEdit()
        self.line_edit_max_y_range = QLineEdit()

        self.layout_y_range = QHBoxLayout()
        self.layout_y_range.addWidget(self.label_y_range)
        self.layout_y_range.addWidget(self.line_edit_min_y_range)
        self.layout_y_range.addWidget(self.line_edit_max_y_range)

        # Reset
        self.button_reset = QPushButton()
        self.button_reset.setText("Reset")

        # Cursor position
        self.line_edit_cursor_pos = QLabel()
        # self.line_edit_cursor_pos.setReadOnly(True)
        self.line_edit_cursor_pos.setText("Pos: 0, 0")

        # Line labels
        self._linelist_window = None
        self._line_labels = []

        self.horizontal_layout.addWidget(self.line_edit_cursor_pos)
        self.horizontal_layout.addStretch()
        # self.horizontal_layout.addLayout(self.layout_x_range)
        # self.horizontal_layout.addLayout(self.layout_y_range)
        self.horizontal_layout.addWidget(self.button_reset)

        self.vertical_layout.addLayout(self.horizontal_layout)

        self.main_widget = QWidget(self)
        self.main_widget.setLayout(self.vertical_layout)

        self.setCentralWidget(self.main_widget)


class PlotSubWindow(UiPlotSubWindow):
    """
    Sub window object responsible for displaying and interacting with plots.
    """
    def __init__(self, *args, **kwargs):
        super(PlotSubWindow, self).__init__(*args, **kwargs)
        self._plots = []
        self._dynamic_axis = None
        self._plot_widget = None
        self._plot_item = None
        self._plots_units = None
        self._rois = []
        self._measure_rois = []
        self._centroid_roi = None
        self._is_selected = True

        DispatchHandle.setup(self)

        self._dynamic_axis = DynamicAxisItem(orientation='top')
        self._plot_widget = pg.PlotWidget(
            axisItems={'top': self._dynamic_axis})
        # self.setCentralWidget(self._plot_widget)
        self.vertical_layout.insertWidget(0, self._plot_widget)

        self._plot_item = self._plot_widget.getPlotItem()
        self._plot_item.showAxis('top', True)
        # Add grids to the plot
        self._plot_item.showGrid(True, True)

        self._setup_connections()

    def _setup_connections(self):
        # Connect cursor position to UI element
        # proxy = pg.SignalProxy(self._plot_item.scene().sigMouseMoved,
        #                        rateLimit=30, slot=self.cursor_moved)
        self._plot_item.scene().sigMouseMoved.connect(self.cursor_moved)
        self.button_reset.clicked.connect(self._reset_view)

    def cursor_moved(self, evt):
        pos = evt

        # Data range
        # flux = self._containers[0].data.value
        # disp = self._containers[0].dispersion.value

        # Plot range
        vb = self._plot_item.getViewBox()

        if self._plot_item.sceneBoundingRect().contains(pos):
            mouse_point = vb.mapSceneToView(pos)
            index = int(mouse_point.x())

            if index > 0 and index < vb.viewRange()[0][1]:
                self.line_edit_cursor_pos.setText("Pos: {0:4.4g}, "
                                                  "{1:4.4g}".format(
                    mouse_point.x(), mouse_point.y())
                    # flux[index], disp[index])
                )

    def _reset_view(self):
        view_box = self._plot_item.getViewBox()
        view_box.autoRange()

    def get_roi_mask(self, layer=None, container=None, roi=None):
        if layer is not None:
            container = self.get_plot(layer)

        if container is None:
            return

        mask_holder = []
        rois = [roi] if roi is not None else self._rois

        for roi in rois:
            # roi_shape = roi.parentBounds()
            # x1, y1, x2, y2 = roi_shape.getCoords()
            x1, x2 = roi.getRegion()

            layer_mask = np.copy(layer._mask)
            mask = (container.dispersion.value >= x1) & \
                   (container.dispersion.value <= x2)
            layer_mask[layer_mask==True] = mask
            mask_holder.append(layer_mask)

        if len(mask_holder) == 0:
            mask_holder.append(container.layer._mask)

        # mask = np.logical_not(reduce(np.logical_or, mask_holder))
        mask = reduce(np.logical_or, mask_holder)
        return mask

    def add_roi(self, *args, **kwargs):
        view_range = self._plot_item.viewRange()
        x_len = (view_range[0][1] - view_range[0][0]) * 0.5
        x_pos = x_len * 0.5 + view_range[0][0]

        def remove():
            self._plot_item.removeItem(roi)
            self._rois.remove(roi)

        roi = LinearRegionItem(values=[x_pos, x_pos + x_len])
        self._rois.append(roi)
        self._plot_item.addItem(roi)

        # Connect the remove functionality
        roi.sigRemoveRequested.connect(remove)

        # Connect events
        Dispatch.on_updated_rois.emit(rois=self._rois)
        roi.sigRemoveRequested.connect(
            lambda: Dispatch.on_updated_rois.emit(rois=self._rois))
        roi.sigRegionChangeFinished.connect(
            lambda: Dispatch.on_updated_rois.emit(rois=self._rois))

    def get_plot(self, layer):
        for plot in self._plots:
            if plot.layer == layer:
                return plot

    def get_all_layers(self):
        return [plot.layer for plot in self._plots]

    def change_units(self, x=None, y=None, z=None):
        for cntr in self._plots:
            cntr.change_units(x, y, z)

        self.set_labels(x_label=x, y_label=y)
        self._plot_item.enableAutoRange()
        self._plot_units = [x, y, z]

    def set_labels(self, x_label='', y_label=''):
        self._plot_item.setLabels(
            left="Flux [{}]".format(
                y_label or str(self._plots[0].layer.units[1])),
            bottom="Wavelength [{}]".format(
                x_label or str(self._plots[0].layer.units[0])))

    def set_visibility(self, layer, show, override=False):
        plot = self.get_plot(layer)

        if plot._visibility_state != [show, show, False]:
            plot.set_visibility(show, show, inactive=False, override=override)

    def update_axis(self, layer=None, mode=None, **kwargs):
        self._dynamic_axis.update_axis(layer, mode, **kwargs)
        self._plot_widget.update()

    def update_plot_item(self):
        self._plot_item.update()

    @DispatchHandle.register_listener("on_update_model")
    def update_plot(self, layer=None, plot=None):
        if layer is not None:
            plot = self.get_plot(layer)

        plot.update()

    def closeEvent(self, event):
        DispatchHandle.tear_down(self)
        super(PlotSubWindow, self).closeEvent(event)

    @DispatchHandle.register_listener("on_add_layer")
    def add_plot(self, layer, window=None):
        if window is not None and window != self:
            logging.warning("Attempted to add container to plot, but sub "
                            "windows do not match.")
            return

        new_plot = LinePlot.from_layer(layer)

        if len(self._plots) == 0:
            self.change_units(new_plot.layer.units[0],
                              new_plot.layer.units[1])
        else:
            new_plot.change_units(*self._plot_units)

        if new_plot.error is not None:
            self._plot_item.addItem(new_plot.error)

        self._plots.append(new_plot)
        self._plot_item.addItem(new_plot.plot)

        self.set_active_plot(new_plot.layer)

        # Make sure the dynamic axis object has access to a layer
        self._dynamic_axis._layer = self._plots[0].layer

        Dispatch.on_added_plot.emit(plot=new_plot, window=window)

    @DispatchHandle.register_listener("on_removed_layer")
    def remove_plot(self, layer, window=None):
        if window is not None and window != self:
            return

        for plot in self._plots:
            if plot.layer == layer:
                self._plot_item.removeItem(plot.plot)

                if plot.error is not None:
                    self._plot_item.removeItem(plot.error)

                self._plots.remove(plot)

    @DispatchHandle.register_listener("on_selected_plot")
    def set_active_plot(self, layer):
        for plot in self._plots:
            if plot.layer == layer:
                plot.set_visibility(True, True, inactive=False)
            else:
                plot.set_visibility(True, False, inactive=True)


#--------  Line lists and line labels handling.

    # The separation of tasks among these methods and the signals
    # that drive them is so far unclear. The 'request linelists'
    # and the 'add_linelists' operations are now in practice
    # fused into a single swoop. It should be more efficient
    # to have the line list ingestion done at constructor time,
    # while only the actual plot should be commanded by the Draw
    # button. This is why the two methods: _request_linelists, and
    # _add_linelists, do not talk directly to each other, but via
    # signals. It's just in preparation for a fix.

    @DispatchHandle.register_listener("on_request_linelists")
    def _request_linelists(self, *args, **kwargs):

        # Find the wavelength range spanned by the spectrum
        # (or spectra) at hand. The range will be used to bracket
        # the set of lines actually read from the line list table(s).

        # increasing dispersion values!
        amin = sys.float_info.max
        amax = 0.0
        for container in self._plots:
            amin = min(amin, container.dispersion.value[0])
            amax = max(amax, container.dispersion.value[-1])

        amin = Quantity(amin, self._plot_units[0])
        amax = Quantity(amax, self._plot_units[0])

        linelist = LineList.ingest(amin, amax)

        Dispatch.on_add_linelists.emit(linelist=linelist)

    @DispatchHandle.register_listener("on_add_linelists")
    def _add_linelists(self, linelist):

        # This is plotting all markers at a fixed height in the
        # screen coordinate system. Still TBD how to do this in
        # the generic case. Maybe derive heights from curve data
        # instead? Make the markers follow the curve ups and downs?
        #
        # Ideally we would like to have the marker's X coordinate
        # pinned down to the plot surface in data value, and the Y
        # coordinate pinned down in screen value. This would make
        # the markers to stay at the same height in the window even
        # when the plot is zoomed. This kind of functionality doesn't
        # seem to be possible under pyqtgraph though. This requires
        # more investigation.

        plot_item = self._plot_item

        # curve = plot_item.curves[0]

        data_range = plot_item.vb.viewRange()
        ymin = data_range[1][0]
        ymax = data_range[1][1]
        height = (ymax - ymin) * 0.75 + ymin

        # column names are defined in the YAML files.
        wave_column = linelist.columns['wavelength']
        id_column = linelist.columns['id']

        for i in range(len(wave_column)):
            marker = LineIDMarker(id_column[i], plot_item, orientation='vertical')

            marker.setPos(wave_column[i], height)

            plot_item.addItem(marker)
            self._line_labels.append(marker)

        plot_item.update()

    @DispatchHandle.register_listener("on_erase_linelabels")
    def erase_linelabels(self, *args, **kwargs):
        for marker in self._line_labels:
            self._plot_item.removeItem(marker)
        self._plot_item.update()


    # The 3 handlers below, and their associated signals, implement
    # the logic that defines the show/hide/dismiss behavior of the
    # line list window. It remains to be seen if it is what users
    # actually want.

    @DispatchHandle.register_listener("on_activated_window")
    def _set_selection_state(self, window):
        self._is_selected = window == self

        if self._linelist_window:
            if self._is_selected:
                self._linelist_window.show()
            else:
                self._linelist_window.hide()

    @DispatchHandle.register_listener("on_show_linelists_window")
    def _show_linelists_window(self, *args, **kwargs):
        if self._is_selected:
            if self._linelist_window is None:
                self._linelist_window = LineListsWindow(str(self))
            self._linelist_window.show()

    @DispatchHandle.register_listener("on_dismiss_linelists_window")
    def _dismiss_linelists_window(self, *args, **kwargs):
        if self._is_selected and self._linelist_window:
            self._linelist_window.hide()
            self._linelist_window = None

