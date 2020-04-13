from astropy.nddata import NDIOMixin


__all__ = ['SpectrumList']


class SpectrumList(list, NDIOMixin):
    """
    A list that is used to hold a list of `~specutils.Spectrum1D` objects

    The primary purpose of this class is to allow loaders to return a list of
    spectra that have different shapes.  For spectra that have the same shape
    but different spectral axes, see `~specutils.SpectrumCollection`.  For
    a spectrum or spectra that all share the same spectral axis, use
    `~specutils.Spectrum1D`.  For more on this topic, see
    :ref:`specutils-representation-overview`.
    """
