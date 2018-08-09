from weakref import WeakKeyDictionary


class LazyProperty:
    def __init__(self, deferred_computation, attr_key=None, 
                 cache_class=WeakKeyDictionary):
        self.deferred_computation = deferred_computation
        self.cache = cache_class()
        self.attr_key = attr_key

    def __get__(self, instance, type=None):
        if instance is None:
            return self

        key = (instance if self.attr_key is None 
               else getattr(instance, self.attr_key))

        if key in self.cache:
            return self.cache[key]
        else:
            result = self.deferred_computation(instance)
            self.cache[key] = result

            return result

    @staticmethod
    def decorate(**kwargs):
        def decorator(func):
            return LazyProperty(func, **kwargs)

        return decorator