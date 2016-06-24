from specviz.interfaces.decorators import data_loader

from astropy.io import ascii

from specviz.core import linelist
from specviz.core.linelist import LineList


def linelist_reader(filename, filter, **kwargs):
    ref = loader_registry.get(filter)

    names_list = []
    start_list = []
    end_list = []
    units_list = []

    for k in range(len((ref.columns))):
        name = ref.columns[k][linelist.COLUMN_NAME]
        names_list.append(name)

        start = ref.columns[k][linelist.COLUMN_START]
        end = ref.columns[k][linelist.COLUMN_END]
        start_list.append(start)
        end_list.append(end)

        if linelist.UNITS_COLUMN in ref.columns[k]:
            units = ref.columns[k][linelist.UNITS_COLUMN]
        else:
            units = ''
        units_list.append(units)

    tab = ascii.read(filename, format = ref.format,
                     names = names_list,
                     col_starts = start_list,
                     col_ends = end_list)

    for k, colname in enumerate(tab.columns):
        tab[colname].unit = units_list[k]

    return LineList(tab)


def linelist_identify(origin, *args, **kwargs):
    """Check whether given filename is a line list.
    """
    return (isinstance(args[0], str) and
            args[0].lower().split('.')[-1] in ['txt', 'dat'])