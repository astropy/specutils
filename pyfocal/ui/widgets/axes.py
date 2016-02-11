import pyqtgraph as pg
import numpy as np
import astropy.constants as const
import astropy.units as u


class DynamicAxisItem(pg.AxisItem):
    """
    """
    def __init__(self, *args, **kwargs):
        super(DynamicAxisItem, self).__init__(*args, **kwargs)
        print("initialized")
        self.mode = 'velocity'
        self._layer = None
        self._redshift = 0.0
        self._ref_wave = 0.0

    def update_axis(self, layer, mode, **kwargs):
        if mode == 'redshift':
            self._redshift = kwargs['redshift']
        elif mode == 'velocity':
            self._ref_wave = kwargs['ref_wave']
        elif mode == 'channel':
            pass

        self._layer = layer or self._layer
        self.mode = mode or self.mode

    def tickStrings(self, values, scale, spacing):
        spatial_unit = self._layer.dispersion.unit

        if self.mode == 'redshift':
            self.setLabel('Redshifted Wavelength [{}]'.format(spatial_unit))

            return [v/(1 + self._redshift)*scale for v in values]

        elif self.mode == 'velocity':
            self.setLabel("Velocity [km/s]", None, None)

            c = const.c.to('{}/s'.format(spatial_unit))
            waves = u.Quantity(np.array(values), spatial_unit)
            ref_wave = u.Quantity(self._ref_wave, spatial_unit)
            v = (waves - ref_wave) / waves * c

            return v.to('km/s').value