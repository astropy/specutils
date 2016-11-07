import os
from collections import OrderedDict

from glue.core import Subset
from glue.viewers.common.qt.data_viewer import DataViewer
from glue.core import message as msg
from glue.utils import nonpartial
from glue.viewers.common.qt.toolbar import BasicToolbar

from specviz.ui.viewer import Viewer
from specviz.core import dispatch
from specviz.core import DispatchHandle

from .viewer_options import OptionsWidget
from .layer_widget import LayerWidget


__all__ = ['SpecVizViewer']


class BaseVizViewer(DataViewer):
    def __init__(self, session, parent=None):
        super(BaseVizViewer, self).__init__(session, parent=parent)

        # Connect the dataview to the specviz messaging system
        DispatchHandle.setup(self)

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
        self._options_widget.set_data(self._layer_widget.layer)
        self._refresh_data()
        return True

    def add_subset(self, subset):
        if subset not in self._layer_widget:
            self._layer_widget.add_layer(subset)
        self._layer_widget.layer = subset
        self._options_widget.set_data(self._layer_widget.layer)
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


class SpecVizViewer(BaseVizViewer):
    LABEL = "SpecViz Viewer"

    def __init__(self, session, parent=None):
        super(SpecVizViewer, self).__init__(session, parent=None)
        # We keep a cache of the specviz data objects that correspond to a given
        # filename - although this could take up a lot of memory if there are
        # many spectra, so maybe this isn't needed
        self._specviz_data_cache = OrderedDict()

        # We set up the specviz viewer and controller as done for the standalone
        # specviz application
        self.viewer = Viewer(hide_plugins=False)
        self.setCentralWidget(self.viewer.main_window)

    def initialize_toolbar(self):
        pass

    def open_data(self, data):
        dispatch.on_add_data.emit(data)

    def _refresh_data(self):
        if self._options_widget.file_att is None:
            print("returning because of options widget")
            return

        if self._layer_widget.layer is None:
            print("returning because of layer widget")
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
        dispatch.on_remove_all_data.emit()

        if not component.categorical:
            return

        filenames = component.labels
        path = '/'.join(component._load_log.path.split('/')[:-1])

        if mask is not None:
            filenames = filenames[mask]

        for filename in filenames:

            if filename in self._specviz_data_cache:
                data = self._specviz_data_cache[filename]
                dispatch.on_add_data.emit(data=data)

            else:
                file_name = str(filename)
                file_path = os.path.join(path, file_name)
                dispatch.on_file_read.emit(file_name=file_path,
                                           file_filter='MOS')

    @DispatchHandle.register_listener('on_added_data')
    def _added_data(self, data):
        print("Adding data")
        filename = data.name
        self._specviz_data_cache[filename] = data
