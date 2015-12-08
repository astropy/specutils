import yaml


class Registry(object):
    """
    Maintains a set of objects.
    """
    def __init__(self):
        self._members = []


class CustomLoaderRegistry(Registry):
    def __init__(self):
        super(CustomLoaderRegistry, self).__init__()
        import os

        for file_name in os.listdir(os.path.join(__file__, 'default_loaders')):
            custom_loader = yaml.load(file_name)
            self._members.append(custom_loader)

    def get(self, filter):
        return [x for x in self._members if x.filter == filter][0]


class YAMLLoader(yaml.YAMLObject):
    yaml_tag = u'!CustomLoader'

    def __init__(self, extension, name, data, uncertainty, mask, meta):
        self.name = name
        self.extension = extension
        self.data = data
        self.uncertainty = uncertainty
        self.mask = mask
        self.meta = meta
        self.filter = "{} (*.{})".format(self.name, self.extension)


loader_factory = CustomLoaderRegistry()