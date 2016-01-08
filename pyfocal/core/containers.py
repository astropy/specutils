from .events import EventHook
from astropy.units import Unit, Quantity


class PlotContainer(object):
    def __init__(self, layer, visible=True, style='line', pen=None):
        self.layer = layer
        self.visible = visible
        self.style = style
        self.pen = pen

        self.plot = None

        self.on_unit_change = EventHook()
        self.on_visibility_change = EventHook()
        self.on_pen_change = EventHook()

    def change_unit(self, x, y=None, z=None):
        self.layer.layer_units = (x, y or self.layer.layer_units[1],
                                  z or self.layer.layer_units[2])

    def change_visible(self, visible_state):
        self.visible = visible_state

    def change_pen(self, new_pen):
        self.pen = new_pen

    @property
    def data(self):
        return self.layer.data

    @property
    def dispersion(self):
        return self.layer.dispersion

    @property
    def units(self):
        return self.layer.layer_units
