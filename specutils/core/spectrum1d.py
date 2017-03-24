from astropy.nddata import NDDataRef


__all__ = ['Spectrum1D']

class Spectrum1D(NDDataRef):
    def __init__(self, *args, **kwargs):
        super(Spectrum1D, self).__init__(*args, **kwargs)

