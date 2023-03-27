"""Defines extension that is used by ASDF for recognizing specutils types."""

__all__ = []


def get_extensions():
    from asdf.extension import ManifestExtension
    from specutils.io.asdf.converters import SpectralAxisConverter, Spectrum1DConverter, SpectrumListConverter

    SPECUTILS_TRANSFORM_CONVERTERS = [SpectralAxisConverter(), Spectrum1DConverter(), SpectrumListConverter()]

    # The order here is important; asdf will prefer to use extensions
    # that occur earlier in the list.
    TRANSFORM_MANIFEST_URIS = [
        "asdf://astropy.org/specutils/manifests/specutils-1.0.0"]

    return [
        ManifestExtension.from_uri(
            uri,
            legacy_class_names=["specutils.io.asdf.extension.SpecutilsExtension"],
            converters=SPECUTILS_TRANSFORM_CONVERTERS)
        for uri in TRANSFORM_MANIFEST_URIS]


def get_resource_mappings():
    from pathlib import Path
    from asdf.resource import DirectoryResourceMapping

    resources_root = Path(__file__).resolve().parent
    if not resources_root.is_dir():
        raise RuntimeError(f"Missing resources directory: {resources_root}")

    return [
        DirectoryResourceMapping(
            resources_root / "schemas",
            "http://astropy.org/schemas/specutils/spectra/"),
        DirectoryResourceMapping(
            resources_root / "manifests",
            "asdf://astropy.org/specutils/manifests/")]
