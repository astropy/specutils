import os

def get_package_data():
    paths = ['coveragerc',
             os.path.join('data', '*fits')]
             os.path.join('data', '*pck')]

    return {'specutils.tests': paths}
