from ..spectra import SpectralRegion


__all__ = ['computation_wrapper']


def computation_wrapper(func, spectrum, region):

    # No region, therefore whole spectrum.
    if region is None:
        return func(spectrum)

    # Single region
    elif isinstance(region, SpectralRegion):
        return func(spectrum, region=region)

    # List of regions
    elif isinstance(region, list):
        return [func(spectrum, region=reg) for reg in region]
