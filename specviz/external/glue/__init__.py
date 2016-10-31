def setup():
    from .data_viewer import SpecVizViewer, MOSVizViewer
    from glue.config import qt_client
    qt_client.add(SpecVizViewer)
    qt_client.add(MOSVizViewer)
