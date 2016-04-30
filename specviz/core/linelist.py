from astropy.table import Table


class LineList(Table):

    @classmethod
    def read(cls, *args, **kwargs):
        from ..interfaces.registries import io_registry
        return io_registry.read(cls, *args, **kwargs)

