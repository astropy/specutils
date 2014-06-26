def get_package_data():
    # Installs the testing data files. Unable to get package_data
    # to deal with a directory hierarchy of files, so just explicitly list.
    return {'specutils.io.tests': ['files/lwp11854.mxlo', 'files/swp02283.mxlo',
                                   'files/*.fits', 'files/*.dat',
                                   'files/*.txt']}
