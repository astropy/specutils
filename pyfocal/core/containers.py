from .events import EventHook
from astropy.units import Unit, Quantity


class PlotContainer(object):
    def __init__(self, layer, unit=None, visible=True, style='line', pen=None):
        self.layer = layer
        self.unit = unit or Unit("")
        self.visible = visible
        self.style = style
        self.pen = pen

        self.plot = None

        self.on_unit_change = EventHook()
        self.on_visibility_change = EventHook()
        self.on_pen_change = EventHook()

    def change_unit(self, new_unit):
        self.unit = Unit(new_unit)

    def change_visible(self, visible_state):
        self.visible = visible_state

    def change_pen(self, new_pen):
        self.pen = new_pen

    @property
    def data(self):
        return Quantity(self.layer.data).to(self.unit).value
