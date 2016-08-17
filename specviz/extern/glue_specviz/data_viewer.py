import os
from collections import OrderedDict

from glue.core import Subset
from glue.viewers.common.qt.data_viewer import DataViewer
from glue.core import message as msg
from glue.utils import nonpartial

from specviz.ui.viewer import Viewer
from specviz.core import Dispatch as SVDispatch
from specviz.core import DispatchHandle as SVDispatchHandle
from mosviz.core import Dispatch as MVDispatch
from mosviz.core import DispatchHandle as MVDispatchHandle

from mosviz.app import App

from .viewer_options import OptionsWidget
from .layer_widget import LayerWidget


__all__ = ['MOSVizViewer']


class BaseVizViewer(DataViewer):
    def __init__(self, session, parent=None):
        super(BaseVizViewer, self).__init__(session, parent=parent)

        # Connect the dataview to the specviz messaging system
        SVDispatchHandle.setup(self)

        # We now set up the options widget. This controls for example which
        # attribute should be used to indicate the filenames of the spectra.
        self._options_widget = OptionsWidget(data_viewer=self)

        # The layer widget is used to select which data or subset to show.
        # We don't use the default layer list, because in this case we want to
        # make sure that only one dataset or subset can be selected at any one
        # time.
        self._layer_widget = LayerWidget()

        # Make sure we update the viewer if either the selected layer or the
        # column specifying the filename is changed.
        self._layer_widget.ui.combo_active_layer.currentIndexChanged.connect(
            nonpartial(self._update_options))
        self._layer_widget.ui.combo_active_layer.currentIndexChanged.connect(
            nonpartial(self._refresh_data))
        self._options_widget.ui.combo_file_attribute.currentIndexChanged.connect(
            nonpartial(self._refresh_data))

    # The following two methods are required by glue - they are used to specify
    # which widgets to put in the bottom left and middle left panel.

    def options_widget(self):
        return self._options_widget

    def layer_view(self):
        return self._layer_widget

    # The following method is required by glue - it is used to subscribe the
    # viewer to various messages sent by glue.

    def register_to_hub(self, hub):

        super(BaseVizViewer, self).register_to_hub(hub)

        hub.subscribe(self, msg.SubsetCreateMessage,
                      handler=self._add_subset)

        hub.subscribe(self, msg.SubsetUpdateMessage,
                      handler=self._update_subset)

        hub.subscribe(self, msg.SubsetDeleteMessage,
                      handler=self._remove_subset)

        hub.subscribe(self, msg.DataUpdateMessage,
                      handler=self._update_data)

    # The following two methods are required by glue - they are what gets called
    # when a dataset or subset gets dragged and dropped onto the viewer.

    def add_data(self, data):
        if data not in self._layer_widget:
            self._layer_widget.add_layer(data)
        self._layer_widget.layer = data
        self._refresh_data()
        return True

    def add_subset(self, subset):
        if subset not in self._layer_widget:
            self._layer_widget.add_layer(subset)
        self._layer_widget.layer = subset
        self._refresh_data()
        return True

    # The following four methods are used to receive various messages related
    # to updates to data or subsets.

    def _update_data(self, message):
        self._refresh_data()

    def _add_subset(self, message):
        self.add_subset(message.subset)

    def _update_subset(self, message):
        self._refresh_data()

    def _remove_subset(self, message):
        if message.subset in self._layer_widget:
            self._layer_widget.remove_layer(message.subset)
        self._refresh_data()

    # When the selected layer is changed, we need to update the combo box with
    # the attributes from which the filename attribute can be selected. The
    # following method gets called in this case.

    def _update_options(self):
        self._options_widget.set_data(self._layer_widget.layer)

    def _refresh_data(self):
        raise NotImplementedError()


class SpecvizViewer(BaseVizViewer):

    LABEL = "SpecViz viewer"

    def __init__(self, session, parent=None):
        super(SpecvizViewer, self).__init__(session, parent=None)
        # We keep a cache of the specviz data objects that correspond to a given
        # filename - although this could take up a lot of memory if there are
        # many spectra, so maybe this isn't needed
        self._specviz_data_cache = OrderedDict()

        # We set up the specviz viewer and controller as done for the standalone
        # specviz application
        self.viewer = Viewer(hide_plugins=True)
        self.setCentralWidget(self.viewer.main_window)

    def _refresh_data(self):
        if self._options_widget.file_att is None:
            return

        if self._layer_widget.layer is None:
            return

        if isinstance(self._layer_widget.layer, Subset):
            subset = self._layer_widget.layer
            cid = subset.data.id[self._options_widget.file_att[0]]
            mask = subset.to_mask(None)
            component = subset.data.get_component(cid)
        else:
            cid = self._layer_widget.layer.id[self._options_widget.file_att[0]]
            mask = None
            component = self._layer_widget.layer.get_component(cid)

        # Clear current data objects in SpecViz
        SVDispatch.on_remove_all_data.emit()

        if not component.categorical:
            return

        filenames = component.labels
        path = '/'.join(component._load_log.path.split('/')[:-1])

        if mask is not None:
            filenames = filenames[mask]

        for filename in filenames:

            if filename in self._specviz_data_cache:
                data = self._specviz_data_cache[filename]
                SVDispatch.on_add_data.emit(data=data)

            else:
                file_name = str(filename)
                file_path = os.path.join(path, file_name)
                SVDispatch.on_file_read.emit(file_name=file_path,
                                           file_filter='MOS')

    @SVDispatchHandle.register_listener('on_added_data')
    def _added_data(self, data):
        filename = data.name
        self._specviz_data_cache[filename] = data


class MOSVizViewer(BaseVizViewer):
    LABEL = "MOSViz viewer"

    def __init__(self, session, parent=None):
        super(MOSVizViewer, self).__init__(session, parent=None)
        # We keep a cache of the mosviz data objects that correspond to a given
        # filename - although this could take up a lot of memory if there are
        # many spectra, so maybe this isn't needed
        self._mosviz_data_cache = OrderedDict()

        self.app = App(full_ui=False)
        self.setCentralWidget(self.app.main_window)

    # def _remove_subset(self, message):
    #     if message.subset in self._layer_widget:
    #         self._layer_widget.remove_layer(message.subset)
    #         del self._mosviz_data_cache[message.subset]
    #
    #     self._refresh_data()

    def _refresh_data(self):
        print(len(self._layer_widget._layers))
        for layer in self._layer_widget._layers:
            if isinstance(layer, Subset):
                mask = layer.to_mask(None)
                data = layer.data
            else:
                mask = None
                data = layer

            columns = []
            col_names = data.components

            for att in col_names:
                cid = data.id[att]
                component = data.get_component(cid)

                if component.categorical:
                    col_data = component.labels[mask]
                else:
                    col_data = component.data[mask]

                if mask is None:
                    col_data = col_data[0]

                if str(att) in ['spectrum1d', 'spectrum2d', 'cutout']:
                    path = '/'.join(component._load_log.path.split('/')[:-1])
                    columns.append([os.path.join(path, x) for x in col_data])
                else:
                    columns.append(col_data)

            catalog_list = []

            for row in zip(*columns):
                catalog_dict = dict(zip([str(x) for x in col_names], row))
                catalog_list.append(catalog_dict)

            if len(catalog_list) == 0:
                continue

            if layer in self._mosviz_data_cache:
                if self._mosviz_data_cache[layer] == catalog_list:
                    continue

            self._mosviz_data_cache[layer] = catalog_list

            print("-"*20, len(self._mosviz_data_cache))

            # Clear current data objects in MOSViz
            MVDispatch.on_remove_all_data.emit()

            for k, v in self._mosviz_data_cache.items():
                MVDispatch.on_file_read.emit(file_name=v,
                                             file_filter='mos-glue')
