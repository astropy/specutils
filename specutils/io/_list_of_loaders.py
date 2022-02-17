import io
from specutils import Spectrum1D
import specutils.io.default_loaders  # noqa
"""
The purpose of this file is to receive a list of loaders from
specutils.spectrum1d.read.list_formats(), format that list
into something that can be used by `automodapi`, and then
set it as the __doc__.
"""


def _list_of_loaders():
    # Receive list of loaders
    list_of_loaders = io.StringIO()
    Spectrum1D.read.list_formats(list_of_loaders)

    # Use the second line (which uses "-" to split the
    # first row from the rest) to create the "=" signs
    # which are used to create the table in .rst files
    split_list = list_of_loaders.getvalue().split("\n")
    line_of_equals = split_list[1].replace("-", "=")
    split_list[1] = line_of_equals

    # Combine elements to create formatted table
    # which can be displayed using `automodapi`
    formatted_table = line_of_equals + "\n"
    for line in split_list:
        formatted_table += line + "\n"
    formatted_table += line_of_equals + "\n"
    return formatted_table


__doc__ = _list_of_loaders()
