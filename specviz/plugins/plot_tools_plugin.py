"""
Manage plot attributes
"""
import os
import logging
from collections import OrderedDict

from astropy.units import Unit

from ..core.comms import dispatch, DispatchHandle
from ..ui.widgets.utils import ICON_PATH
from ..ui.widgets.plugin import Plugin
from ..ui.widgets.dialogs import TopAxisDialog, UnitChangeDialog


class PlotToolsPlugin(Plugin):
    """
    UI plugin to manage plot attributes of the various layers
    """
    name = "Plot Tools"
    location = "hidden"
    _all_categories = {}

    def setup_ui(self):
        self._top_axis_dialog = TopAxisDialog()
        self._unit_change_dialog = UnitChangeDialog()

        # Add an roi
        self.button_add_roi = self.add_tool_bar_actions(
            name="ROI",
            description='Add ROI',
            icon_path=os.path.join(ICON_PATH, "Merge Vertical-48.png"),
            category=('Selections', 4),
            priority=1,
            callback=dispatch.on_add_roi.emit,
            enabled=False)

        # Change top axis
        self.button_axis_change = self.add_tool_bar_actions(
            name="Top Axis",
            description='Change top axis',
            icon_path=os.path.join(ICON_PATH, "Globe Earth-48.png"),
            category=('Options', 2),
            callback=self._top_axis_dialog.exec_,
            enabled=False)

        # Change plot units
        self.button_unit_change = self.add_tool_bar_actions(
            name="Units",
            description='Change plot units',
            icon_path=os.path.join(ICON_PATH, "Generic Text-48.png"),
            category='Options',
            priority=1,
            callback=self._show_unit_change_dialog,
            enabled=False)

        self.button_line_labels = self.add_tool_bar_actions(
            name="Line Labels",
            description='Add line labels',
            icon_path=os.path.join(ICON_PATH, "Label-48.png"),
            category='Selections',
            callback=dispatch.on_show_linelists_window.emit,
            enabled=False)

        self.button_plot_settings = self.add_tool_bar_actions(
            name="Plot Settings",
            description='Edit visual plot settings',
            icon_path=os.path.join(ICON_PATH, "Settings-50.png"),
            category='Options',
            callback=dispatch.on_show_linelists_window.emit,
            enabled=False,
            menu=OrderedDict([
                ('Plot style', OrderedDict([
                    ('Line', lambda: self._set_plot_style(mode='line')),
                    ('Scatter', lambda: self._set_plot_style(mode='scatter')),
                    ('Histogram', lambda: self._set_plot_style(mode='histogram'))
                ])),
                ('Line width', OrderedDict([
                    ('1', lambda: self._set_plot_style(line_width=1)),
                    ('2', lambda: self._set_plot_style(line_width=2)),
                    ('3', lambda: self._set_plot_style(line_width=3))
                ])),
                ('Show Errors', ['checkable', lambda x: self._toggle_errors(x)])
            ])
        )

    def setup_connections(self):
        # On accept, change the displayed axis
        self._top_axis_dialog.accepted.connect(
            self._update_axis)

    def _show_unit_change_dialog(self):
        # Populate the text fields with the current units
        self._unit_change_dialog.line_edit_flux_unit.setText("{}".format(self.current_layer.unit))
        self._unit_change_dialog.line_edit_disp_unit.setText("{}".format(self.current_layer.dispersion_unit))

        if self._unit_change_dialog.exec_():
            x_text = self._unit_change_dialog.disp_unit
            y_text = self._unit_change_dialog.flux_unit

            x_unit = y_unit = None

            try:
                x_unit = Unit(x_text) if x_text else None
            except ValueError as e:
                logging.error(e)

            try:
                y_unit = Unit(y_text) if y_text else None
            except ValueError as e:
                logging.error(e)

            self.active_window.change_units(x_unit, y_unit)
            self.active_window.update_plot_item()

    def _update_axis(self):
        if self.active_window is None:
            return

        if len(self.active_window._plots) > 0:
            layer = self.active_window._plots[0].layer

            self.active_window.update_axis(
                layer,
                self._top_axis_dialog.combo_box_axis_mode.currentIndex(),
                redshift=self._top_axis_dialog.redshift,
                ref_wave=self._top_axis_dialog.ref_wave)
        else:
            logging.warning("Active window does not have any plots.")

    def _set_plot_style(self, **kwargs):
        if self.active_window is not None:
            self.active_window.set_plot_style(self.current_layer, **kwargs)

    def _toggle_errors(self, state):
        if self.active_window is not None:
            layer = self.current_layer
            current_window = self.active_window
            current_window.disable_errors = not state
            current_window.set_active_plot(layer)

    @DispatchHandle.register_listener("on_activated_window")
    def toggle_enabled(self, window):
        if window:
            self.button_axis_change.setEnabled(True)
            self.button_unit_change.setEnabled(True)
            self.button_add_roi.setEnabled(True)
            # self.button_line_labels.setEnabled(True)
            self.button_plot_settings.setEnabled(True)
        else:
            self.button_axis_change.setEnabled(False)
            self.button_unit_change.setEnabled(False)
            self.button_add_roi.setEnabled(False)
            # self.button_line_labels.setEnabled(False)
            self.button_plot_settings.setEnabled(False)
