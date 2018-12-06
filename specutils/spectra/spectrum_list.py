from astropy.nddata import NDIOMixin


__all__ = ['SpectrumList']


class SpectrumList(list, NDIOMixin):
    """
    A list that is used to hold a list of Spectrum1D objects

    The primary purpose of this class is to allow loaders to return a list of
    heterogenous spectra that do have a spectral axis of the same length.
    """
