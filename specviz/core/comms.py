"""
Object Event Handling

The singleton `Dispatch` object manages the
set of `EventNode` events. Handlers or **listeners** are attached
to `EventNode`s. The `DispatchHandle` decorator is used to
decorate classes that handle events.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import inspect
import traceback
import logging
from functools import wraps


class EventNode(object):
    """An event

    An event is defined by the arguments the listeners
    of the event expect to be given.

    Parameters
    ----------
    args: [arg, ...]
        The list of keyword arguments that the event
        will provide to its listeners
    """
    def __init__(self, *args):
        self._args = args
        self.__handlers = []

    def __iadd__(self, other):
        """
        Add a listener via the addition operator
        """
        self.__handlers.append(other)
        return self

    def __isub__(self, other):
        """
        Remove a listener via the subtraction operator
        """
        self.__handlers.remove(other)
        return self

    def emit(self, *args, **kwargs):
        """
        Call the hanlders of this event

        Parameters
        ----------
        args: [arg, ...]
            The keyword arguments being provided to the event.

        kwargs: {arg: value, ...}
            The keyword/value pairs passed to the listeners.

        Raises
        ------
        ValueError
            An keyword is being passed that does not belong
            to this event.
        """
        if len(args) != len(self._args) and not set(kwargs.keys()).issubset(
                set(self._args)):
            raise ValueError("Unknown keyword in event emit arguments.")

        for handler in self.__handlers:
            # logging.info("Sending message from: '{}'".format(handler))
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
    def register_event(self, name, args=None):
        """
        Add an `EventNode` to the list of possible events

        Parameters
        ----------
        name: str
            The name of the event.

        args: [arg, ...]
            The list of keyword arguments this event will pass
            to its handlers.
        """
        args = args or []

        if not hasattr(self, name):
            setattr(self, name, EventNode(*args))
        else:
            logging.warning("Event '{}' already exists. Please use a "
                            "different name.".format(name))

    def register_listener(self, name, func):
        """
        Add a listener to an event

        Parameters
        ----------
        name: str
            The event name to add the listener to.

        func: function
            The function that will be called when the is emitted
        """
        if hasattr(self, name):
            call_func = getattr(self, name)
            call_func += func
        else:
            logging.warning("No such event: {}. Event must be registered "
                            "before listeners can be assigned.".format(name))

    def unregister_listener(self, name, func):
        """
        Remove a listener from an event

        Parameters
        ----------
        name: str
            The event from wich the listener should be removed.

        func: function
            The function to be removed
        """
        if hasattr(self, name):
            call_func = getattr(self, name)
            call_func -= func
        else:
            logging.warning("No such event: {}.".format(name))


class DispatchHandle(object):
    """
    Interface for allowing classes to use decorators to define event
    listeners. Otherwise, classes would have to define all listeners in the
    `init` function using

    >>> dispatch.register_listener("<event_name>", <class_method>)
    """
    @staticmethod
    def setup(inst):
        """
        Register all methods decorated by `register_listener`
        """
        logging.info("Dispatch is now watching: {}".format(inst))
        members = inspect.getmembers(inst, predicate=inspect.ismethod)

        for func_name, func in members:
            if hasattr(func, 'wrapped'):
                if func.wrapped:
                    for name in func.event_names:
                        dispatch.register_listener(name, func)

    @staticmethod
    def tear_down(inst):
        """
        Remove all registered methods from their events
        """
        logging.info("Dispatch has stopped watching: {}".format(inst))
        members = inspect.getmembers(inst, predicate=inspect.ismethod)

        for func_name, func in members:
            if hasattr(func, 'wrapped'):
                if func.wrapped:
                    for name in func.event_names:
                        dispatch.unregister_listener(name, func)

    @staticmethod
    def register_listener(*args):
        """
        Decorate event listeners
        """
        def decorator(func):
            func.wrapped = True
            func.event_names = args

            @wraps(func)
            def wrapper(*args, **kwargs):
                try:
                    return func(*args, **kwargs)
                except:
                    logging.error(
                        "Exception in '{}':\n{}".format(func.__name__,
                                                        traceback.format_exc())
                    )
            return wrapper
        return decorator


# Register application-wide events
dispatch = Dispatch()
dispatch.register_event("on_activated_window", args=["window"])

dispatch.register_event("on_added_data", args=["data"])
dispatch.register_event("on_added_window", args=["layer", "window"])
dispatch.register_event("on_added_plot", args=["plot", "window"])
dispatch.register_event("on_added_layer", args=["layer"])
dispatch.register_event("on_added_to_window", args=["layer", "window"])

dispatch.register_event("on_show_linelists_window")
dispatch.register_event("on_dismiss_linelists_window")
dispatch.register_event("on_request_linelists")
dispatch.register_event("on_plot_linelists", args=["table_views"])
dispatch.register_event("on_erase_linelabels")

dispatch.register_event("on_removed_data", args=["data"])
dispatch.register_event("on_removed_plot", args=["layer", "window"])
dispatch.register_event("on_removed_layer", args=["layer", "window"])
dispatch.register_event("on_removed_model", args=["model", "layer"])
dispatch.register_event("on_removed_from_window", args=["layer", "window"])

dispatch.register_event("on_updated_layer", args=["layer"])
dispatch.register_event("on_updated_model", args=["model"])
dispatch.register_event("on_updated_plot", args=["plot", "layer"])
dispatch.register_event("on_updated_rois", args=["rois"])
dispatch.register_event("on_updated_stats", args=["stats", "layer"])

dispatch.register_event("on_selected_plot", args=["layer", "checked_state"])
dispatch.register_event("on_selected_window", args=["window"])
dispatch.register_event("on_selected_layer", args=["layer_item"])
dispatch.register_event("on_selected_model", args=["model_item"])

dispatch.register_event("on_clicked_layer", args=["layer_item"])
dispatch.register_event("on_changed_layer", args=["layer_item"])
dispatch.register_event("on_changed_model", args=["model_item"])

dispatch.register_event("on_add_data", args=["data"])
dispatch.register_event("on_add_model", args=["layer"])
dispatch.register_event("on_add_window", args=["data", "window"])
dispatch.register_event("on_add_layer", args=["window", "layer", "from_roi"])
dispatch.register_event("on_add_roi", args=[])

dispatch.register_event("on_update_model", args=["layer"])

dispatch.register_event("on_remove_data", args=["data"])
dispatch.register_event("on_remove_layer", args=["layer"])
dispatch.register_event("on_remove_model", args=["model"])
dispatch.register_event("on_remove_all_data")

dispatch.register_event("on_file_open", args=["file_name"])
dispatch.register_event("on_file_read", args=["file_name", "file_filter"])

dispatch.register_event("on_status_message", args=["message", "timeout"])
