import os

def get_package_data():
    paths = ['coveragerc',
             os.path.join('data', '*fits')]

    return {'specutils.tests': paths}
