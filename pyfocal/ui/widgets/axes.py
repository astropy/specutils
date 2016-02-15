import pyqtgraph as pg
import numpy as np
import astropy.constants as const
import astropy.units as u


class DynamicAxisItem(pg.AxisItem):
    """
    """
    def __init__(self, *args, **kwargs):
        super(DynamicAxisItem, self).__init__(*args, **kwargs)
        self.mode = 0
        self.supported_modes = ['velocity', 'redshift', 'channel']
        self._layer = None
        self.redshift = 0.0
        self.ref_wave = 0.0

    def update_axis(self, layer, mode, **kwargs):
        if mode == 1:
            self.redshift = kwargs['redshift']
        elif mode == 0:
            self.ref_wave = kwargs['ref_wave']
        elif mode == 2:
            pass

        self._layer = layer or self._layer
        self.mode = self.supported_modes[mode]

    def generateDrawSpecs(self, p):
        if self.mode == 'channel':
            mn, mx = self.range[0], self.range[1]
            self.setLabel("Channel", None, None)
            data = self._xdata[(self._xdata > mn) & (self._xdata < mx)]
            self.setTicks([
                [(v, str(i)) for i, v in enumerate(data)][::len(
                    data)/10]
            ])
        return super(DynamicAxisItem, self).generateDrawSpecs(p)

    def tickStrings(self, values, scale, spacing):
        if self._layer is None:
            return super(DynamicAxisItem, self).tickStrings(values, scale,
                                                            spacing)

        spatial_unit = self._layer.dispersion.unit

        if self.mode == 'redshift':
            self.setLabel('Redshifted Wavelength [{}]'.format(spatial_unit))

            return [v / (1 + self.redshift) * scale for v in values]

        elif self.mode == 'velocity':
            self.setLabel("Velocity [km/s]", None, None)

            c = const.c.to('{}/s'.format(spatial_unit))

            waves = u.Quantity(np.array(values), spatial_unit)
            waves[np.isnan(waves)] = 0.0

            ref_wave = u.Quantity(self.ref_wave, spatial_unit)

            v = (waves - ref_wave) / waves * c
            v = v.to('km/s').value
            v[np.isnan(v)] = 0.0

            return list(v)