def setup():
    from .data_viewer import SpecvizViewer, MOSVizViewer
    from glue.config import qt_client
    qt_client.add(SpecvizViewer)
    qt_client.add(MOSVizViewer)
