from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import inspect
import traceback
import logging
from functools import wraps


class EventNode(object):
    def __init__(self, *args):
        self._args = args
        self.__handlers = []

    def __iadd__(self, other):
        self.__handlers.append(other)
        return self

    def __isub__(self, other):
        self.__handlers.remove(other)
        return self

    def emit(self, *args, **kwargs):
        if len(args) != len(self._args) and not set(kwargs.keys()).issubset(
                set(self._args)):
            raise ValueError("Unknown keyword in event emit arguments.")

        for handler in self.__handlers:
            if hasattr(handler, 'self'):
                handler(handler.self, *args, **kwargs)
            else:
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
    def register_event(cls, name, args=None):
        args = args or []

        if not hasattr(cls, name):
            setattr(cls, name, EventNode(*args))
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

    @classmethod
    def unregister_listener(cls, name, func):
        if hasattr(cls, name):
            call_func = getattr(cls, name)
            call_func -= func
        else:
            logging.warning("No such event: {}.".format(name))


class DispatchHandle(object):
    """
    Interface for allowing classes to use decorators to define event
    listeners. Otherwise, classes would have to define all listeners in the
    `init` function using

    >>> Dispatch.register_listener("<event_name>", <class_method>)
    """
    @staticmethod
    def setup(inst):
        logging.info("Dispatch is now watching: {}".format(inst))
        members = inspect.getmembers(inst, predicate=inspect.ismethod)

        for func_name, func in members:
            if hasattr(func, 'wrapped'):
                if func.wrapped:
                    for name in func.event_names:
                        Dispatch.register_listener(name, func)

    @staticmethod
    def tear_down(inst):
        logging.info("Dispatch has stopped watching: {}".format(inst))
        members = inspect.getmembers(inst, predicate=inspect.ismethod)

        for func_name, func in members:
            if hasattr(func, 'wrapped'):
                if func.wrapped:
                    for name in func.event_names:
                        Dispatch.unregister_listener(name, func)

    @staticmethod
    def register_listener(*args):
        def decorator(func):
            func.wrapped = True
            func.event_names = args

            @wraps(func)
            def wrapper(*args, **kwargs):
                try:
                    return func(*args, **kwargs)
                except:
                    logging.error(traceback.format_exc())
            return wrapper
        return decorator


Dispatch.register_event("on_added_data", args=["data"])
Dispatch.register_event("on_added_window", args=["data", "window"])
Dispatch.register_event("on_added_plot", args=["container", "window"])
Dispatch.register_event("on_added_model", args=["model", "layer"])
Dispatch.register_event("on_added_layer", args=["layer"])
Dispatch.register_event("on_added_to_window", args=["layer", "window"])

Dispatch.register_event("on_removed_data", args=["data"])
Dispatch.register_event("on_removed_plot", args=["layer", "window"])
Dispatch.register_event("on_removed_layer", args=["layer"])
Dispatch.register_event("on_removed_model", args=["model", "layer"])
Dispatch.register_event("on_removed_from_window", args=["layer", "window"])

Dispatch.register_event("on_updated_layer", args=["layer"])
Dispatch.register_event("on_updated_model", args=["model"])
Dispatch.register_event("on_updated_plot", args=["container", "layer"])
Dispatch.register_event("on_updated_roi", args=["roi", "measured_rois"])
Dispatch.register_event("on_updated_stats", args=["stats", "layer"])

Dispatch.register_event("on_selected_plot", args=["layer"])
Dispatch.register_event("on_selected_window", args=["window"])
Dispatch.register_event("on_selected_layer", args=["layer_item"])
Dispatch.register_event("on_selected_model", args=["model_item"])

Dispatch.register_event("on_clicked_layer", args=["layer_item"])
Dispatch.register_event("on_changed_model", args=["model_item"])
