import copy
import warnings

import astropy.units as u
from astropy.nddata import StdDevUncertainty
from astropy.utils.exceptions import AstropyUserWarning

from ...spectra import Spectrum, SpectrumList
from ..parsing_utils import read_fileobj_or_asdftree
from ..registers import data_loader

# Optional Roman GDPS package
# if the user has this installed, roman spectral files serialize into pydantic datamodels
try:
    import roman_gdps
except ImportError:
    roman_gdps = None


def _identify_roman_1d(*args, mode: str, units: str=None, **kwargs) -> bool:
    """Identify a Roman SSC spectral file

    Parameters
    ----------
    mode : str
        the spectral mode
    units : str, optional
        the flux units substring to check, by default None

    Returns
    -------
    bool
        whether the input matches the given conditions
    """
    with read_fileobj_or_asdftree(*args, **kwargs) as af:
        # fail if not a roman asdf file
        if "asdf_library" not in af.tree or "roman" not in af.tree:
            return False

        # normalize to standard python dicts
        if roman_gdps:
            roman = af.tree["roman"].model_dump()
        else:
            roman = af.tree.get("roman", {})

        # set base condition
        meta = roman.get("meta", {})
        base = (isinstance(roman, dict) and
                (meta.get('optical_element') or meta.get('instr_optical_element')) in ("GRISM", "PRISM") and
                "data" in roman and isinstance(roman["data"], dict)
                )

        # check specific conditions
        if mode == 'single':
            return base and "flux" in roman["data"]
        if mode == 'multi':
            return base and "flux" not in roman["data"] and "1d_spectra" in roman["meta"].get("title", "") and units in roman['meta']['unit_flux']


def identify_1d_single_source(origin, *args, **kwargs):
    """Check if input is a Roman single source 1d extracted spectrum"""

    return _identify_roman_1d(*args, mode='single', **kwargs)

def identify_1d_combined(origin, *args, **kwargs):
    """Check if input is a set of Roman 1d combined extracted spectra"""

    return _identify_roman_1d(*args, mode='multi', units='W m**(-2)', **kwargs)


def identify_1d_individual(origin, *args, **kwargs):
    """Check if input is a set of Roman 1d individual extracted spectra"""

    return _identify_roman_1d(*args, mode='multi', units='DN/s', **kwargs)


def _load_roman_spectrum(roman: dict, source: str) -> Spectrum:
    """Load a single Roman spectrum

    Create and return a Spectrum object from a Roman SSC input dictionary.

    Parameters
    ----------
    roman : dict
        A collection of Roman spectral data and metadata.
    source : str
        A source identifier string.

    Returns
    -------
    Spectrum
        the Roman Spectrum object
    """

    meta = copy.deepcopy(roman["meta"])
    meta['source_id'] = source
    data = roman["data"][source] if source else roman["data"]
    flux = data['flux'] * u.Unit(meta["unit_flux"])
    flux_err  = StdDevUncertainty(data['flux_error'])
    wavelength = data['wl'] * u.Unit(meta["unit_wl"])
    return Spectrum(spectral_axis=wavelength, flux=flux, uncertainty=flux_err, meta=meta)


def _load_roman_first(file_obj, source=None, **kwargs):
    """Load the Roman spectrum for the first source found"""
    with read_fileobj_or_asdftree(file_obj, **kwargs) as af:
        roman = af["roman"]
        if source is None:
            source = next(iter(roman['data']))
            warnings.warn(f"source not specified. Loading first source found: {source}.",
                        AstropyUserWarning)
        elif source not in roman['data']:
            raise ValueError(f"Invalid source! Source {source} is not spectra.")

        return _load_roman_spectrum(roman, source)


def _lazy_load_roman(file_obj, **kwargs):
    """Lazy loader for SpectrumList"""

    # read the file and get the sources
    with read_fileobj_or_asdftree(file_obj, **kwargs) as af:
        roman = af["roman"]
        sources = list(roman["data"].keys())
        source_idx_map = dict(zip(sources, range(len(sources))))

    # define the lazy loader
    def _loader(i: int) -> Spectrum:
        source = sources[i]
        with read_fileobj_or_asdftree(file_obj, **kwargs) as af2:
            roman2 = af2["roman"]
            return _load_roman_spectrum(roman2, source)

    # create the lazy object
    cache_size = kwargs.get("cache_size", None)
    sl = SpectrumList.from_lazy(
        length=len(sources), loader=_loader, cache_size=cache_size, labels=sources
    )

    # set the sources ids as alternate indices
    sl.set_id_map(source_idx_map)

    return sl


def _load_roman_multisource(file_obj, **kwargs) -> SpectrumList:
    """Load all Roman spectra into a SpectrumList"""

    spectra = SpectrumList()
    with read_fileobj_or_asdftree(file_obj, **kwargs) as af:
        roman = af["roman"]
        meta = roman["meta"]
        sources = list(roman['data'].keys())

        # set the alternate ids to roman source ids
        source_idx_map = dict(zip(sources, range(len(sources))))
        meta['source_map'] = source_idx_map
        spectra.set_id_map(source_idx_map)
        # load the spectra
        for source in roman["data"]:
            spectrum = _load_roman_spectrum(roman, source)
            spectra.append(spectrum)

    return spectra


@data_loader(
    "Roman 1d spectra", identifier=identify_1d_single_source, dtype=Spectrum,
    extensions=['asdf'], priority=10,
)
def roman_1d_spectrum_loader(file_obj, **kwargs):
    """Load a Roman single source 1d extracted spectrum

    Loader for Roman 1d extracted spectral data in ASDF format.

    Parameters
    ----------
    file_obj: str, file-like, or HDUList
          FITS file name, object (provided from name by Astropy I/O Registry),
          or HDUList (as resulting from astropy.io.fits.open()).

    Returns
    -------
    Spectrum
        The spectrum contained in the file.
    """
    with read_fileobj_or_asdftree(file_obj, **kwargs) as af:
        roman = af["roman"]
        meta = roman["meta"]
        data = roman["data"]
        flux = data['flux'] * u.Unit(meta["unit_flux"])
        flux_err  = StdDevUncertainty(data['flux_error'])
        wavelength = data['wl'] * u.Unit(meta["unit_wl"])
        return Spectrum(spectral_axis=wavelength, flux=flux, uncertainty=flux_err, meta=meta)



@data_loader(
    "Roman 1d combined",
    identifier=identify_1d_combined,
    dtype=Spectrum,
    extensions=["asdf"],
    priority=10,
    autogenerate_spectrumlist=False,
)
def roman_1d_combined(file_obj, source: str = None, **kwargs):
    """Load Roman 1d combined extracted spectra

    Loads a set of extracted spectra from a Roman 1d combined ASDF file.
    If not source is specified, it loads the first source found in the data dict.

    Parameters
    ----------
    file_obj
        the input file or file object
    source : str, optional
        the source id string, by default None

    Returns
    -------
    Spectrum
        The Roman spectrum object
    """
    return _load_roman_first(file_obj, source=source, **kwargs)


@data_loader(
    "Roman 1d combined",
    identifier=identify_1d_combined,
    dtype=SpectrumList,
    extensions=["asdf"],
    priority=10,
    force=True,
    lazy_loader=_lazy_load_roman,
)
def roman_1d_combined_list(file_obj, **kwargs):
    """Load all Roman 1d combined extracted spectra

    Loads a set of extracted spectra from a Roman 1d combined ASDF file,
    as a list of Spectrum objects.

    Parameters
    ----------
    file_obj
        the input file or file object

    Returns
    -------
    SpectrumList
        A set of Roman spectrum object
    """
    return _load_roman_multisource(file_obj, **kwargs)


@data_loader(
    "Roman 1d individual",
    identifier=identify_1d_individual,
    dtype=Spectrum,
    extensions=["asdf"],
    priority=10,
    autogenerate_spectrumlist=False,
)
def roman_1d_individual(file_obj, source: str=None, **kwargs):
    """Load Roman 1d individual extracted spectra

    Loads a set of extracted spectra from a Roman 1d individual ASDF file.
    If not source is specified, it loads the first source found in the data dict.

    Parameters
    ----------
    file_obj
        the input file or file object
    source : str, optional
        the source id string, by default None

    Returns
    -------
    Spectrum
        The Roman spectrum object
    """
    return _load_roman_first(file_obj, source=source, **kwargs)


@data_loader(
    "Roman 1d individual",
    identifier=identify_1d_individual,
    dtype=SpectrumList,
    extensions=["asdf"],
    priority=10,
    force=True,
    lazy_loader=_lazy_load_roman,
)
def roman_1d_individual_list(file_obj, **kwargs):
    """Load all Roman 1d individual extracted spectra

    Loads a set of extracted spectra from a Roman 1d individual ASDF file,
    as a list of Spectrum objects.

    Parameters
    ----------
    file_obj
        the input file or file object
    source : str, optional
        the source id string, by default None

    Returns
    -------
    SpectrumList
        A set of Roman spectrum object
    """
    return _load_roman_multisource(file_obj, **kwargs)


