from setuptools import setup


setup(
    name='SpecViz',
    version='0.1a1',
    packages=['specviz', 'specviz.ui', 'specviz.ui.qt', 'specviz.ui.widgets',
              'specviz.ui.widgets.plots', 'specviz.core', 'specviz.analysis',
              'specviz.interfaces',
              'specviz.third_party',
              'specviz.third_party.py_expression_eval',
              'specviz.third_party.qtpy'],
    url='http://specviz.readthedocs.org',
    license='',
    author='Nicholas Earl, Ivo Busko, Pey Lian Lim',
    author_email='nearl@stsci.edu',
    description='An interactive astronomical analysis tool.',
    include_package_data=True,
    install_requires=[
        'astropy',
        'numpy>=1.10',
        'PyYAML',
        'pyqtgraph',
        'scipy'
    ],
    entry_points={
        'console_scripts': [
            'specviz = specviz.app:main'
        ]
    }
)
