from importlib.util import find_spec
import pkg_resources

entry_points = []
for entry_point in pkg_resources.iter_entry_points('pytest11'):
    entry_points.append(entry_point.name)

if "asdf_schema_tester" not in entry_points and find_spec('asdf') is not None:
	pytest_plugins = ['asdf.tests.schema_tester']

