from ..spectra import SpectralRegion


__all__ = ['computation_wrapper']


def computation_wrapper(func, spectrum, region, **kwargs):

    # No region, therefore whole spectrum.
    if region is None:
        return func(spectrum, **kwargs)

    # Single region
    elif isinstance(region, SpectralRegion):
        return func(spectrum, region=region, **kwargs)

    # List of regions
    elif isinstance(region, list):
        return [func(spectrum, region=reg, **kwargs) for reg in region]
