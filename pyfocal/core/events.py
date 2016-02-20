from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging
from functools import wraps


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
        if args:
            raise ValueError("Undefined arguments in event emit.")

        if not set(kwargs.keys()).issubset(set(self._kwargs.keys())):
            raise ValueError("Unknown keyword in event emit arguments.")

        for handler in self.__handlers:
            if hasattr(handler, 'self'):
                handler(handler.self, **kwargs)
            else:
                handler(**kwargs)

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
    def register_event(cls, name, kwargs=()):
        if not hasattr(cls, name):
            setattr(cls, name, EventHook(**{k: None for k in kwargs}))
        else:
            logging.warning("Event '{}' already exists. Please use a "
                            "different name.".format(name))

    @classmethod
    def register_listener(cls, name, func):
        if hasattr(cls, name):
            call_func = getattr(cls, name)
            call_func += func
        else:
            logging.warning("No such event: {}. Event must be registered "
                            "before listeners can be assigned.".format(name))


class DispatchHandle(object):
    """
    Interface for allowing classes to use decorators to define event
    listeners. Otherwise, classes would have to define all listeners in the
    `init` function using

    >>> Dispatch.register_listener("<event_name>", <class_method>)
    """
    @staticmethod
    def setup(inst):
        for func_name in dir(inst):
            func = getattr(inst, func_name)

            if hasattr(func, 'wrapped'):
                func().self = inst

                if func.wrapped:
                    Dispatch.register_listener(func.event_name, func())

    @staticmethod
    def register_listener(name):
        def decorator(func):
            func.wrapped = True
            func.event_name = name
            print(func.__class__)

            @wraps(func)
            def wrapper(*args, **kwargs):
                return func
            return wrapper
        return decorator


Dispatch.register_event("on_add_data", kwargs=("data",))
Dispatch.register_event("on_update_layer", kwargs=("layer",))
Dispatch.register_event("on_update_model", kwargs=("data",))
Dispatch.register_event("on_update_stats", kwargs=("data",))
Dispatch.register_event("on_set_plot_active", kwargs=("data",))


