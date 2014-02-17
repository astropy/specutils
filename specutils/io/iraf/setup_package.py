def get_package_data():
    # Installs the testing data files. Unable to get package_data
    # to deal with a directory hierarchy of files, so just explicitly list.
    return {'specutils.io.iraf.tests': ['data/idexample']}
