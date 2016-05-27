from ..ui.widgets.plugin import Plugin
from ..third_party.qtpy.QtWidgets import *
from ..third_party.qtpy.QtCore import *
from ..core.comms import Dispatch, DispatchHandle
from ..interfaces.managers import layer_manager
from ..analysis import statistics

import logging


class StatisticsPlugin(Plugin):
    name = "Statistics"

    def __init__(self, *args, **kwargs):
        super(StatisticsPlugin, self).__init__(*args, **kwargs)
        self._current_window = None
        self._current_layer_item = None

    def setup_ui(self):
        self.layout_vertical.setContentsMargins(11, 11, 11, 11)

        # Setup form layout
        self.layout_form = QFormLayout()
        self.layout_form.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)
        self.layout_form.setFormAlignment(Qt.AlignJustify | Qt.AlignTop)
        self.layout_form.setContentsMargins(1, 1, 1, 12)
        self.layout_form.setSpacing(6)

        self.layout_vertical.addLayout(self.layout_form)

        # Setup labels
        self.label_current_layer = QLabel(self)
        self.label_current_layer.setText("Current Layer")

        self.line_edit_current_layer = QLineEdit(self)
        sizePolicy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.line_edit_current_layer.sizePolicy().hasHeightForWidth())
        self.line_edit_current_layer.setSizePolicy(sizePolicy)
        self.line_edit_current_layer.setStyleSheet("QLineEdit{background: #DDDDDD;}")
        self.line_edit_current_layer.setReadOnly(True)

        self.layout_form.setWidget(0, QFormLayout.LabelRole,
                                   self.label_current_layer)
        self.layout_form.setWidget(0, QFormLayout.FieldRole,
                                   self.line_edit_current_layer)

        # Setup tabs
        self.tab_widget_stats = QTabWidget(self)
        self.layout_vertical.addWidget(self.tab_widget_stats)

        # Setup basic tab
        self.tab_basic = QWidget()
        self.tab_widget_stats.addTab(self.tab_basic, "Basic")

        self.layout_vertical_tab_basic = QVBoxLayout(self.tab_basic)
        self.layout_vertical_tab_basic.setContentsMargins(11, 11, 11, 11)
        self.layout_vertical_tab_basic.setSpacing(6)

        self.layout_form_tab_basic = QFormLayout()
        self.layout_form_tab_basic.setFieldGrowthPolicy(
            QFormLayout.ExpandingFieldsGrow)
        self.layout_form_tab_basic.setFormAlignment(Qt.AlignRight |
                                                    Qt.AlignTop |
                                                    Qt.AlignTrailing)
        self.layout_form_tab_basic.setContentsMargins(1, 1, 1, 1)
        self.layout_form_tab_basic.setSpacing(6)

        self.layout_vertical_tab_basic.addLayout(self.layout_form_tab_basic)

        # Setup basic tab labels
        self.label_mean = QLabel(self.tab_basic)
        self.label_mean.setText("Mean")

        self.line_edit_mean = QLineEdit(self.tab_basic)
        self.line_edit_mean.setStyleSheet(
            "QLineEdit{background: #DDDDDD;}")
        self.line_edit_mean.setReadOnly(True)
        self.layout_form_tab_basic.setWidget(0, QFormLayout.LabelRole,
                                             self.label_mean)
        self.layout_form_tab_basic.setWidget(0, QFormLayout.FieldRole,
                                             self.line_edit_mean)

        self.label_median = QLabel(self.tab_basic)
        self.label_median.setText("Median")

        self.line_edit_median = QLineEdit(self.tab_basic)
        self.line_edit_median.setStyleSheet(
            "QLineEdit{background: #DDDDDD;}")
        self.line_edit_median.setReadOnly(True)
        self.layout_form_tab_basic.setWidget(1, QFormLayout.LabelRole,
                                             self.label_median)
        self.layout_form_tab_basic.setWidget(1, QFormLayout.FieldRole,
                                             self.line_edit_median)

        self.label_std_dev = QLabel(self.tab_basic)
        self.label_std_dev.setText("Standard Deviation")

        self.line_edit_std_dev = QLineEdit(self.tab_basic)
        self.line_edit_std_dev.setStyleSheet(
            "QLineEdit{background: #DDDDDD;}")
        self.line_edit_std_dev.setReadOnly(True)
        self.layout_form_tab_basic.setWidget(2, QFormLayout.LabelRole,
                                             self.label_std_dev)
        self.layout_form_tab_basic.setWidget(2, QFormLayout.FieldRole,
                                             self.line_edit_std_dev)

        self.label_total = QLabel(self.tab_basic)
        self.label_total.setText("Total")

        self.line_edit_total = QLineEdit(self.tab_basic)
        self.line_edit_total.setStyleSheet(
            "QLineEdit{background: #DDDDDD;}")
        self.line_edit_total.setReadOnly(True)
        self.layout_form_tab_basic.setWidget(3, QFormLayout.LabelRole,
                                             self.label_total)
        self.layout_form_tab_basic.setWidget(3, QFormLayout.FieldRole,
                                             self.line_edit_total)

        self.label_data_point_count = QLabel(self.tab_basic)
        self.label_data_point_count.setText("Data Point Count")

        self.line_edit_data_point_count = QLineEdit(self.tab_basic)
        self.line_edit_data_point_count.setStyleSheet(
            "QLineEdit{background: #DDDDDD;}")
        self.line_edit_data_point_count.setReadOnly(True)
        self.layout_form_tab_basic.setWidget(4, QFormLayout.LabelRole,
                                             self.label_data_point_count)
        self.layout_form_tab_basic.setWidget(4, QFormLayout.FieldRole,
                                             self.line_edit_data_point_count)

        # Measured tab setup
        self.tab_measured = QWidget()
        self.tab_widget_stats.addTab(self.tab_measured, "Measured")

        self.layout_vertical_tab_measured = QVBoxLayout(self.tab_measured)
        self.layout_vertical_tab_measured.setContentsMargins(11, 11, 11, 11)
        self.layout_vertical_tab_measured.setSpacing(6)

        self.layout_form_tab_measured = QFormLayout()
        self.layout_form_tab_measured.setFieldGrowthPolicy(
            QFormLayout.ExpandingFieldsGrow)
        self.layout_form_tab_measured.setFormAlignment(Qt.AlignRight |
                                                       Qt.AlignTop |
                                                       Qt.AlignTrailing)
        self.layout_form_tab_measured.setContentsMargins(1, 1, 1, 1)
        self.layout_form_tab_measured.setSpacing(6)

        self.layout_vertical_tab_measured.addLayout(self.layout_form_tab_measured)

        # Measured tab labels
        self.label_equivalent_width = QLabel(self.tab_measured)
        self.label_equivalent_width.setText("Equivalent Width")

        self.line_edit_equivalent_width = QLineEdit(self.tab_measured)
        self.line_edit_equivalent_width.setStyleSheet(
            "QLineEdit{background: #DDDDDD;}")
        self.line_edit_equivalent_width.setReadOnly(True)

        self.layout_form_tab_measured.setWidget(0, QFormLayout.LabelRole,
                                                self.label_equivalent_width)
        self.layout_form_tab_measured.setWidget(0, QFormLayout.FieldRole,
                                                self.line_edit_equivalent_width)

        self.label_centroid = QLabel(self.tab_measured)
        self.label_centroid.setText("Centroid")

        self.line_edit_centroid = QLineEdit(self.tab_measured)
        self.line_edit_centroid.setStyleSheet(
            "QLineEdit{background: #DDDDDD;}")
        self.line_edit_centroid.setReadOnly(True)

        self.layout_form_tab_measured.setWidget(1, QFormLayout.LabelRole,
                                                self.label_centroid)
        self.layout_form_tab_measured.setWidget(1, QFormLayout.FieldRole,
                                                self.line_edit_centroid)

        self.label_flux = QLabel(self.tab_measured)
        self.label_flux.setText("Flux")

        self.line_edit_flux = QLineEdit(self.tab_measured)
        self.line_edit_flux.setStyleSheet(
            "QLineEdit{background: #DDDDDD;}")
        self.line_edit_flux.setReadOnly(True)

        self.layout_form_tab_measured.setWidget(2, QFormLayout.LabelRole,
                                                self.label_flux)
        self.layout_form_tab_measured.setWidget(2, QFormLayout.FieldRole,
                                                self.line_edit_flux)

        self.label_mean_continuum = QLabel(self.tab_measured)
        self.label_mean_continuum.setText("Mean Continuum")

        self.line_edit_continuum = QLineEdit(self.tab_measured)
        self.line_edit_continuum.setStyleSheet(
            "QLineEdit{background: #DDDDDD;}")
        self.line_edit_continuum.setReadOnly(True)

        self.layout_form_tab_measured.setWidget(3, QFormLayout.LabelRole,
                                                self.label_mean_continuum)
        self.layout_form_tab_measured.setWidget(3, QFormLayout.FieldRole,
                                                self.line_edit_continuum)

        # Add warning label
        self.label_measured_error = QLabel()
        self.label_measured_error.setText("You must have at least three ROIs "
                                          "on the plot")
        self.label_measured_error.setWordWrap(True)
        self.label_measured_error.setStyleSheet("""
        QLabel {
            color: #a94442;
            background-color: #f2dede;
            padding: 10px;
            border: 1px solid #ebccd1;
            border-radius: 4px;
        }""")

        self.layout_vertical_tab_measured.addWidget(self.label_measured_error)

        self.layout_vertical.addStretch()

    def setup_connections(self):
        pass

    @DispatchHandle.register_listener("on_selected_window")
    def set_window(self, window=None):
        self._current_window = window

    @DispatchHandle.register_listener("on_selected_layer")
    def set_layer(self, layer_item=None):
        self._current_layer_item = layer_item
        current_layer = self._current_layer_item.data(0, Qt.UserRole)
        self.line_edit_current_layer.setText(current_layer.name)

    @DispatchHandle.register_listener("on_updated_roi")
    def update_basic_stats(self, roi=None, measure_rois=None):
        if roi is None:
            return

        if self._current_window is None or self._current_layer_item is None:
            logging.warning(
                "No window or layer item provided; cannot update statistics.")
            return

        current_layer = self._current_layer_item.data(0, Qt.UserRole)

        # Set the active tab to basic
        self.tab_widget_stats.setCurrentIndex(0)

        mask = self._current_window.get_roi_mask(layer=current_layer)

        if mask is None:
            values = current_layer.data
        else:
            values = current_layer.data[mask[current_layer._mask]]

        stat_dict = statistics.stats(values)

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

    @DispatchHandle.register_listener("on_updated_roi")
    def update_measured_states(self, roi=None, measure_rois=None):
        if measure_rois is None:
            return

        if self._current_window is None or self._current_layer_item is None:
            logging.warning(
                "No window or layer item provided; cannot update statistics.")
            return

        current_layer = self._current_layer_item.data(0, Qt.UserRole)

        # Set the active tab to measured
        self.tab_widget_stats.setCurrentIndex(1)

        cont1_mask = self._current_window.get_roi_mask(layer=current_layer,
                                                       roi=measure_rois[0])
        cont1_data = current_layer.data[cont1_mask]
        cont1_stat_dict = statistics.stats(cont1_data)

        cont2_mask = self._current_window.get_roi_mask(layer=current_layer,
                                                       roi=measure_rois[2])
        cont2_data = current_layer.data[cont2_mask]
        cont2_stat_dict = statistics.stats(cont2_data)

        line_mask = self._current_window.get_roi_mask(layer=current_layer,
                                                      roi=measure_rois[1])

        line = layer_manager.copy(current_layer)

        ew, flux, avg_cont = statistics.eq_width(cont1_stat_dict,
                                                 cont2_stat_dict,
                                                 line,
                                                 mask=line_mask)
        cent = statistics.centroid(line - avg_cont, mask=line_mask)

        stat_dict = {"eq_width": ew, "centroid": cent, "flux": flux,
                     "avg_cont": avg_cont}

        self.line_edit_equivalent_width.setText("{0:4.4g}".format(
            float(stat_dict['eq_width'].value)))
        self.line_edit_centroid.setText("{0:5.5g}".format(
            float(stat_dict['centroid'].value)))
        self.line_edit_flux.setText("{0:4.4g}".format(
            float(stat_dict['flux'].value)))
        self.line_edit_continuum.setText("{0:4.4g}".format(
            float(stat_dict['avg_cont'].value)))