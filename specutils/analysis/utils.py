"""
A module for internal utilities for the analysis sub-package.  Not meant for
public API consumption.
"""

from ..spectra import SpectralRegion


__all__ = ['computation_wrapper']


def computation_wrapper(func, spectrum, region, **kwargs):
    """
    Applies a computation across either a whole spectrum or a bunch of regions.
    """

    # No region, therefore whole spectrum.
    if region is None:
        return func(spectrum, **kwargs)

    # Single region
    elif isinstance(region, SpectralRegion):
        return func(spectrum, regions=region, **kwargs)

    # List of regions
    elif isinstance(region, list):
        return [func(spectrum, regions=reg, **kwargs) for reg in region]
