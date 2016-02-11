from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


class EventHook(object):
    def __init__(self):
        self.__handlers = []

    def __iadd__(self, other):
        self.__handlers.append(other)
        return self

    def __isub__(self, other):
        self.__handlers.remove(other)
        return self

    def emit(self, *args, **kwargs):
        for handler in self.__handlers:
            handler(*args, **kwargs)

    def clear(self):
        """
        Removes all handlers from object.
        """
        self.__handlers = []
