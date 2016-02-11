from setuptools import setup


setup(
    name='Pyfocal',
    version='0.1a1',
    packages=['pyfocal', 'pyfocal.ui', 'pyfocal.ui.qt',
              'pyfocal.ui.qt.uic.source.qdarkstyle', 'pyfocal.ui.widgets',
              'pyfocal.ui.widgets.plots', 'pyfocal.core', 'pyfocal.analysis',
              'pyfocal.interfaces',
              'pyfocal.third_party',
              'pyfocal.third_party.py_expression_eval',
              'pyfocal.third_party.qtpy'],
    url='http://pyfocal.readthedocs.org',
    license='',
    author='Nicholas Earl, Ivo Busko, Pey Lian Lim',
    author_email='nmearl@protonmail.com',
    description='An interactive astronomical analysis tool.',
    include_package_data=True,
    install_requires=[
        'astropy>=1.1',
        'numpy>=1.10',
        'PyYAML',
        'pyqtgraph',
        'scipy'
    ],
    entry_points={
        'console_scripts': [
            'pyfocal = pyfocal.app:main'
        ]
    }
)
