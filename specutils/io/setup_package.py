def get_package_data():
    # Installs the testing data files. Unable to get package_data
    # to deal with a directory hierarchy of files, so just explicitly list.
    return {'specutils.io.tests': ['files/lwp11854.mxlo', 'files/swp02283.mxlo',
                                   'files/1d.fits', 'files/AAO.fits', 'files/AAO_11.txt', 'files/TRES.dat',
                                   'files/TRES.fits', 'files/UVES.fits', 'files/lwp11854.mxlo',
                                   'files/swp02283.mxlo']
            }
