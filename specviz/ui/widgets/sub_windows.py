from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging
from functools import reduce

from .axes import DynamicAxisItem
from ...third_party.qtpy.QtWidgets import *
from ...third_party.qtpy.QtGui import *
from ...third_party.qtpy.QtCore import *
from ...core.comms import Dispatch, DispatchHandle
from .region_items import LinearRegionItem

import numpy as np
import pyqtgraph as pg

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
pg.setConfigOptions(antialias=False)


class UiPlotSubWindow(QMainWindow):
    def __init__(self, *args, **kwargs):
        super(UiPlotSubWindow, self).__init__(*args, **kwargs)

        self.vertical_layout = QVBoxLayout()
        self.horizontal_layout = QHBoxLayout()

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
        self.line_edit_cursor_pos.setText("None, None")

        self.horizontal_layout.addWidget(self.line_edit_cursor_pos)
        # self.horizontal_layout.addLayout(self.layout_x_range)
        # self.horizontal_layout.addLayout(self.layout_y_range)
        # self.horizontal_layout.addWidget(self.button_reset)
        self.horizontal_layout.addStretch()

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
        self._containers = []
        self._dynamic_axis = None
        self._plot_widget = None
        self._plot_item = None
        self._plots_units = None
        self._rois = []
        self._measure_rois = []
        self._centroid_roi = None

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
                self.line_edit_cursor_pos.setText("{0:4.4g}, {1:4.4g}".format(
                    mouse_point.x(), mouse_point.y())
                    # flux[index], disp[index])
                )

    def get_roi_mask(self, layer=None, container=None, roi=None):
        if layer is not None:
            container = self.get_container(layer)

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

    def get_container(self, layer):
        for container in self._containers:
            if container.layer == layer:
                return container

    def change_units(self, x=None, y=None, z=None):
        for cntr in self._containers:
            cntr.change_units(x, y, z)

        self.set_labels(x_label=x, y_label=y)
        self._plot_item.enableAutoRange()
        self._plot_units = [x, y, z]

    def set_labels(self, x_label='', y_label=''):
        self._plot_item.setLabels(
            left="Flux [{}]".format(
                y_label or str(self._containers[0].layer.units[1])),
            bottom="Wavelength [{}]".format(
                x_label or str(self._containers[0].layer.units[0])))

    def set_visibility(self, layer, show, override=False):
        for container in self._containers:
            if container.layer == layer:
                if container._visibility_state != [show, show, False]:
                    container.set_visibility(show, show, inactive=False,
                                             override=override)

    def update_axis(self, layer=None, mode=None, **kwargs):
        self._dynamic_axis.update_axis(layer, mode, **kwargs)
        self._plot_widget.update()

    def closeEvent(self, event):
        DispatchHandle.tear_down(self)
        super(PlotSubWindow, self).closeEvent(event)

    @DispatchHandle.register_listener("on_added_plot")
    def add_container(self, container, window):
        if window != self:
            logging.warning("Attempted to add container to plot, but sub "
                            "windows do not match.")
            return

        # User tries to plot before loading file, do nothing
        if not hasattr(container, 'layer'):
            return

        if len(self._containers) == 0:
            self.change_units(container.layer.units[0],
                              container.layer.units[1])
        else:
            container.change_units(*self._plot_units)

        if container.error is not None:
            self._plot_item.addItem(container.error)

        self._containers.append(container)
        self._plot_item.addItem(container.plot)

        self.set_active_plot(container.layer)

        # Make sure the dynamic axis object has access to a layer
        self._dynamic_axis._layer = self._containers[0].layer

    @DispatchHandle.register_listener("on_removed_plot")
    def remove_container(self, layer, window=None):
        if window is not None and window != self:
            return

        for container in [x for x in self._containers]:
            if container.layer == layer:
                self._plot_item.removeItem(container.plot)

                if container.error is not None:
                    self._plot_item.removeItem(container.error)

                self._containers.remove(container)

    @DispatchHandle.register_listener("on_selected_plot")
    def set_active_plot(self, layer):
        for container in self._containers:
            if container.layer == layer:
                container.set_visibility(True, True, inactive=False)
            else:
                container.set_visibility(True, False, inactive=True)

