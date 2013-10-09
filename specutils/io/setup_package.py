def get_package_data():
    # Installs the testing data files. Unable to get package_data
    # to deal with a directory hierarchy of files, so just explicitly list.
    return {'specutils.io.tests': ['t/lwp11854.mxlo', 't/swp02283.mxlo',
                                   't/1d.fits', 't/AAO.fits', 't/AAO_11.txt', 't/TRES.dat',
                                   't/TRES.fits', 't/UVES.fits', 't/lwp11854.mxlo',
                                   't/swp02283.mxlo']
            }
