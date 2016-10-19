from ..ui.widgets.plugin import Plugin
from ..third_party.qtpy.QtWidgets import *
from ..third_party.qtpy.QtCore import *
from ..core.comms import dispatch, DispatchHandle
from ..analysis import statistics
from ..third_party.qtpy.QtGui import *

import logging
import pyqtgraph as pg
import numpy as np
import astropy.units as u
from functools import reduce


LINE_EDIT_CSS = "QLineEdit {background: #DDDDDD; border: 1px solid #cccccc;}"


class StatisticsPlugin(Plugin):
    name = "Statistics"
    location = "left"

    def setup_ui(self):
        UiStatisticsPlugin(self)

    def setup_connections(self):
        pass

    @DispatchHandle.register_listener("on_updated_rois", "on_selected_layer")
    def update_statistics(self, rois=None, *args, **kwargs):
        if rois is None:
            if self.active_window is not None:
                rois = self.active_window._rois
            else:
                rois = []

        current_layer = self._current_layer

        if self.active_window is None or current_layer is None:
            logging.info(
                "No window or layer item provided; cannot update statistics.")

            # Clear statistics information
            for att in self.__dict__:
                if 'line_edit' in att:
                    self.__dict__[att].setText("")

            return

        # Set current layer name text
        self.line_edit_current_layer.setText(current_layer.name)

        mask = self.active_window.get_roi_mask(layer=current_layer)

        if mask is None:
            values = current_layer.data
        else:
            values = np.ma.array(current_layer.data, mask=~mask)

        stat_dict = statistics.stats(values.compressed())

        self.line_edit_mean.setText("{0:4.4g}".format(
            stat_dict['mean'].value))
        self.line_edit_median.setText("{0:4.4g}".format(
            stat_dict['median'].value))
        self.line_edit_std_dev.setText("{0:4.4g}".format(
            stat_dict['stddev'].value))
        self.line_edit_total.setText("{0:4.4g}".format(
            stat_dict['total'].value))
        self.line_edit_data_point_count.setText("{0:4.4g}".format(
            stat_dict['npoints']))

        # Calculate measured statistics if there are three rois
        if len(rois) < 3:
            # So that the rois are not updating all the time, reset the
            # colors of the rois when the number has *just* fallen below 3
            if self.label_measured_error.isHidden():
                [x.setBrush(QColor(0, 0, 255, 50)) for x in rois]
                [x.update() for x in rois]

            self.label_measured_error.show()
            return
        else:
            [x.setBrush(QColor(0, 0, 255, 50)) for x in rois]
            [x.update() for x in rois]

            self.label_measured_error.hide()

        roi_data_sets = []
        roi_masks = []

        for roi in rois:
            mask = self.active_window.get_roi_mask(layer=current_layer,
                                                   roi=roi)
            roi_masks.append(mask)
            values = np.ma.array(current_layer.data, mask=~mask)
            roi_data_sets.append(values)

        # Always make the ROI that's over the greatest absolute data value
        # orange
        # roi_data_sets, rois, roi_masks = zip(*sorted(
        #     zip(roi_data_sets, rois, roi_masks),
        #     key=lambda x: np.max(np.abs(x[0]))))

        rois[-1].setBrush(pg.mkBrush(QColor(255, 69, 0, 50)))
        rois[-1].update()

        cont1_stat_dict = statistics.stats(roi_data_sets[0].compressed().value)
        cont2_stat_dict = statistics.stats(
            np.concatenate([x.compressed().value for x in roi_data_sets[:-1]]))

        line = current_layer

        ew, flux, avg_cont = statistics.eq_width(cont1_stat_dict,
                                                 cont2_stat_dict,
                                                 line,
                                                 mask=roi_masks[-1])

        cent = statistics.centroid(flux=line.data.compressed().value - avg_cont,
                                   wave=line.dispersion.compressed().value,
                                   mask=roi_masks[-1])

        stat_dict = {"eq_width": ew, "centroid": cent, "flux": flux,
                     "avg_cont": avg_cont}

        self.line_edit_equivalent_width.setText("{0:4.4g}".format(
            float(stat_dict['eq_width'])))
        self.line_edit_centroid.setText("{0:5.5g}".format(
            float(stat_dict['centroid'])))
        self.line_edit_flux.setText("{0:4.4g}".format(
            float(stat_dict['flux'])))
        self.line_edit_continuum.setText("{0:4.4g}".format(
            float(stat_dict['avg_cont'])))


class UiStatisticsPlugin:
    def __init__(self, plugin):
        plugin.layout_vertical.setContentsMargins(11, 11, 11, 11)

        # Setup form layout
        plugin.layout_form = QFormLayout()
        plugin.layout_form.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)
        plugin.layout_form.setFormAlignment(Qt.AlignJustify | Qt.AlignTop)
        plugin.layout_form.setContentsMargins(1, 1, 1, 12)
        plugin.layout_form.setSpacing(6)

        plugin.layout_vertical.addLayout(plugin.layout_form)

        # Setup labels
        plugin.label_current_layer = QLabel(plugin)
        plugin.label_current_layer.setText("Current Layer")

        plugin.line_edit_current_layer = QLineEdit(plugin)
        sizePolicy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            plugin.line_edit_current_layer.sizePolicy().hasHeightForWidth())
        plugin.line_edit_current_layer.setSizePolicy(sizePolicy)
        plugin.line_edit_current_layer.setStyleSheet(LINE_EDIT_CSS)
        plugin.line_edit_current_layer.setReadOnly(True)

        plugin.layout_form.setWidget(0, QFormLayout.LabelRole,
                                     plugin.label_current_layer)
        plugin.layout_form.setWidget(0, QFormLayout.FieldRole,
                                     plugin.line_edit_current_layer)

        # Setup tabs
        plugin.tab_widget_stats = QTabWidget(plugin)
        plugin.layout_vertical.addWidget(plugin.tab_widget_stats)

        # Setup basic tab
        plugin.tab_basic = QWidget()
        plugin.tab_widget_stats.addTab(plugin.tab_basic, "Basic")

        plugin.layout_vertical_tab_basic = QVBoxLayout(plugin.tab_basic)
        plugin.layout_vertical_tab_basic.setContentsMargins(11, 11, 11, 11)
        plugin.layout_vertical_tab_basic.setSpacing(6)

        plugin.layout_form_tab_basic = QFormLayout()
        plugin.layout_form_tab_basic.setFieldGrowthPolicy(
            QFormLayout.ExpandingFieldsGrow)
        plugin.layout_form_tab_basic.setFormAlignment(Qt.AlignRight |
                                                      Qt.AlignTop |
                                                      Qt.AlignTrailing)
        plugin.layout_form_tab_basic.setContentsMargins(1, 1, 1, 1)
        plugin.layout_form_tab_basic.setSpacing(6)

        plugin.layout_vertical_tab_basic.addLayout(plugin.layout_form_tab_basic)

        # Setup basic tab labels
        plugin.label_mean = QLabel(plugin.tab_basic)
        plugin.label_mean.setText("Mean")

        plugin.line_edit_mean = QLineEdit(plugin.tab_basic)
        plugin.line_edit_mean.setStyleSheet(LINE_EDIT_CSS)
        plugin.line_edit_mean.setReadOnly(True)
        plugin.layout_form_tab_basic.setWidget(0, QFormLayout.LabelRole,
                                               plugin.label_mean)
        plugin.layout_form_tab_basic.setWidget(0, QFormLayout.FieldRole,
                                               plugin.line_edit_mean)

        plugin.label_median = QLabel(plugin.tab_basic)
        plugin.label_median.setText("Median")

        plugin.line_edit_median = QLineEdit(plugin.tab_basic)
        plugin.line_edit_median.setStyleSheet(
            LINE_EDIT_CSS)
        plugin.line_edit_median.setReadOnly(True)
        plugin.layout_form_tab_basic.setWidget(1, QFormLayout.LabelRole,
                                               plugin.label_median)
        plugin.layout_form_tab_basic.setWidget(1, QFormLayout.FieldRole,
                                               plugin.line_edit_median)

        plugin.label_std_dev = QLabel(plugin.tab_basic)
        plugin.label_std_dev.setText("Std. Dev.")

        plugin.line_edit_std_dev = QLineEdit(plugin.tab_basic)
        plugin.line_edit_std_dev.setStyleSheet(
            LINE_EDIT_CSS)
        plugin.line_edit_std_dev.setReadOnly(True)
        plugin.layout_form_tab_basic.setWidget(2, QFormLayout.LabelRole,
                                               plugin.label_std_dev)
        plugin.layout_form_tab_basic.setWidget(2, QFormLayout.FieldRole,
                                               plugin.line_edit_std_dev)

        plugin.label_total = QLabel(plugin.tab_basic)
        plugin.label_total.setText("Total")

        plugin.line_edit_total = QLineEdit(plugin.tab_basic)
        plugin.line_edit_total.setStyleSheet(
            LINE_EDIT_CSS)
        plugin.line_edit_total.setReadOnly(True)
        plugin.layout_form_tab_basic.setWidget(3, QFormLayout.LabelRole,
                                               plugin.label_total)
        plugin.layout_form_tab_basic.setWidget(3, QFormLayout.FieldRole,
                                               plugin.line_edit_total)

        plugin.label_data_point_count = QLabel(plugin.tab_basic)
        plugin.label_data_point_count.setText("Data Point Count")

        plugin.line_edit_data_point_count = QLineEdit(plugin.tab_basic)
        plugin.line_edit_data_point_count.setStyleSheet(
            LINE_EDIT_CSS)
        plugin.line_edit_data_point_count.setReadOnly(True)
        plugin.layout_form_tab_basic.setWidget(4, QFormLayout.LabelRole,
                                               plugin.label_data_point_count)
        plugin.layout_form_tab_basic.setWidget(4, QFormLayout.FieldRole,
                                               plugin.line_edit_data_point_count)

        # Measured tab setup
        plugin.tab_measured = QWidget()
        plugin.tab_widget_stats.addTab(plugin.tab_measured, "Measured")

        plugin.layout_vertical_tab_measured = QVBoxLayout(plugin.tab_measured)
        plugin.layout_vertical_tab_measured.setContentsMargins(11, 11, 11, 11)
        plugin.layout_vertical_tab_measured.setSpacing(6)

        plugin.layout_form_tab_measured = QFormLayout()
        plugin.layout_form_tab_measured.setFieldGrowthPolicy(
            QFormLayout.ExpandingFieldsGrow)
        plugin.layout_form_tab_measured.setFormAlignment(Qt.AlignRight |
                                                         Qt.AlignTop |
                                                         Qt.AlignTrailing)
        plugin.layout_form_tab_measured.setContentsMargins(1, 1, 1, 1)
        plugin.layout_form_tab_measured.setSpacing(6)

        plugin.layout_vertical_tab_measured.addLayout(plugin.layout_form_tab_measured)

        # Measured tab labels
        plugin.label_equivalent_width = QLabel(plugin.tab_measured)
        plugin.label_equivalent_width.setText("Equivalent Width")

        plugin.line_edit_equivalent_width = QLineEdit(plugin.tab_measured)
        plugin.line_edit_equivalent_width.setStyleSheet(
            LINE_EDIT_CSS)
        plugin.line_edit_equivalent_width.setReadOnly(True)

        plugin.layout_form_tab_measured.setWidget(0, QFormLayout.LabelRole,
                                                  plugin.label_equivalent_width)
        plugin.layout_form_tab_measured.setWidget(0, QFormLayout.FieldRole,
                                                  plugin.line_edit_equivalent_width)

        plugin.label_centroid = QLabel(plugin.tab_measured)
        plugin.label_centroid.setText("Centroid")

        plugin.line_edit_centroid = QLineEdit(plugin.tab_measured)
        plugin.line_edit_centroid.setStyleSheet(
            LINE_EDIT_CSS)
        plugin.line_edit_centroid.setReadOnly(True)

        plugin.layout_form_tab_measured.setWidget(1, QFormLayout.LabelRole,
                                                  plugin.label_centroid)
        plugin.layout_form_tab_measured.setWidget(1, QFormLayout.FieldRole,
                                                  plugin.line_edit_centroid)

        plugin.label_flux = QLabel(plugin.tab_measured)
        plugin.label_flux.setText("Flux")

        plugin.line_edit_flux = QLineEdit(plugin.tab_measured)
        plugin.line_edit_flux.setStyleSheet(
            LINE_EDIT_CSS)
        plugin.line_edit_flux.setReadOnly(True)

        plugin.layout_form_tab_measured.setWidget(2, QFormLayout.LabelRole,
                                                  plugin.label_flux)
        plugin.layout_form_tab_measured.setWidget(2, QFormLayout.FieldRole,
                                                  plugin.line_edit_flux)

        plugin.label_mean_continuum = QLabel(plugin.tab_measured)
        plugin.label_mean_continuum.setText("Mean Continuum")

        plugin.line_edit_continuum = QLineEdit(plugin.tab_measured)
        plugin.line_edit_continuum.setStyleSheet(
            LINE_EDIT_CSS)
        plugin.line_edit_continuum.setReadOnly(True)

        plugin.layout_form_tab_measured.setWidget(3, QFormLayout.LabelRole,
                                                  plugin.label_mean_continuum)
        plugin.layout_form_tab_measured.setWidget(3, QFormLayout.FieldRole,
                                                  plugin.line_edit_continuum)

        # Add warning label
        plugin.label_measured_error = QLabel()
        plugin.label_measured_error.setText("You must have at least three ROIs "
                                          "on the plot")
        plugin.label_measured_error.setWordWrap(True)
        plugin.label_measured_error.setStyleSheet("""
        QLabel {
            color: #a94442;
            background-color: #f2dede;
            padding: 10px;
            border: 1px solid #ebccd1;
            border-radius: 4px;
        }""")

        plugin.layout_vertical_tab_measured.addWidget(plugin.label_measured_error)
        plugin.layout_vertical_tab_measured.addStretch()

        plugin.layout_vertical.addStretch()

        # Set size of plugin. Setting this seems to screw with `QPushButton`
        # visual formatting
        plugin.setMinimumSize(plugin.sizeHint())
