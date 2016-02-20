from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging


class EventHook(object):
    def __init__(self, **kwargs):
        self._kwargs = kwargs
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


class Dispatch(object):
    """
    Central communications object for all events.
    """
    @classmethod
    def register_event(cls, name, args=()):
        if not hasattr(cls, name):
            setattr(cls, name, EventHook())
        else:
            logging.warning("Event '{}' already exists. Please use a "
                            "different name.".format(name))

    @classmethod
    def register_listener(cls, name):
        def decorator(func):
            if hasattr(cls, name):
                call_func = getattr(cls, name)
                call_func += func
            else:
                logging.warning("No such event: {}. Event must be registered "
                                "before listeners can be assigned.".format(name))
            return func
        return decorator


Dispatch.register_event("on_update_data")
Dispatch.register_event("on_update_layer")
Dispatch.register_event("on_update_model_layer")

