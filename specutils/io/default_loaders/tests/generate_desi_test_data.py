"""
Generate small test versions of DESI files::

    python -m specutils.io.default_loaders.tests.generate_desi_test_data

The resulting files try to be self-consistent, *i.e.*, as close as possible
to a genuine DESI file with N spectra, where N is of order 10.
"""
import os
import sys
import warnings
from argparse import ArgumentParser
import numpy as np
from astropy.io import fits
from astropy.utils.data import download_file


desi_url = 'https://data.desi.lbl.gov'
downloads = {'healpix_coadd': 'public/edr/spectro/redux/fuji/healpix/sv3/dark/260/26065/coadd-sv3-dark-26065.fits',
             'healpix_spectra': 'public/edr/spectro/redux/fuji/healpix/sv3/dark/260/26065/spectra-sv3-dark-26065.fits',
             'tile_coadd': 'public/edr/spectro/redux/fuji/tiles/cumulative/169/20210419/coadd-5-169-thru20210419.fits',
             'tile_spectra': 'public/edr/spectro/redux/fuji/tiles/cumulative/169/20210419/spectra-5-169-thru20210419.fits'}


def _options():
    """Parse command-line options.
    """
    prsr = ArgumentParser(description='Generate small test versions of DESI files.')
    prsr.add_argument('-l', '--local-path', metavar='PATH', dest='path',
                      help='Find local files in PATH instead of downloading.')
    prsr.add_argument('-o', '--overwrite', action='store_true',
                      help='Overwrite any existing files.')
    prsr.add_argument('-r', '--rows', default=5, type=int, metavar='N',
                      help='Save N spectra from the original file (default %(default)s).')
    prsr.add_argument('output', metavar='PATH', help='Write files to PATH.')
    return prsr.parse_args()


def main():
    """Entry point for command-line scripts.

    Returns
    -------
    :class:`int`
        An integer suitable for passing to :func:`sys.exit`.
    """
    options = _options()
    files = list()
    for key in downloads:
        if options.path is None:
            filename = download_file(f"{desi_url}/{downloads[key]}", cache=True)
        else:
            filename = f"{options.path}/{downloads[key]}"
        if os.path.isfile(filename):
            files.append(filename)
        else:
            warnings.warn(f"Could not find {filename}!")
    for f in files:
        print(f)
        b = os.path.basename(f)
        with fits.open(f, mode='readonly') as hdulist:
            new_hdulist = hdulist.copy()
            for hdu in new_hdulist:
                try:
                    print(hdu.header['EXTNAME'])
                except KeyError:
                    print("PRIMARY")
                if hdu.header['NAXIS'] == 0:
                    pass
                elif hdu.header['EXTNAME'] == 'FIBERMAP' or hdu.header['EXTNAME'] == 'EXP_FIBERMAP' or hdu.header['EXTNAME'] == 'SCORES':
                    if hdu.header['EXTNAME'] == 'FIBERMAP':
                        targetids = hdu.data['TARGETID'][0:options.rows]
                    if hdu.header['EXTNAME'] == 'EXP_FIBERMAP':
                        exp_index = None
                        for t in targetids:
                            if exp_index is None:
                                exp_index = np.where(hdu.data['TARGETID'] == t)[0]
                            else:
                                exp_index = np.append(exp_index, np.where(hdu.data['TARGETID'] == t)[0])
                        hdu.data = hdu.data[exp_index]
                    else:
                        hdu.data = hdu.data[0:options.rows]
                    # hdu.add_checksum()
                elif hdu.header['NAXIS'] == 1:
                    pass
                elif hdu.header['NAXIS'] == 2:
                    hdu.data = hdu.data[0:options.rows, :]
                    # hdu.add_checksum()
                elif hdu.header['NAXIS'] == 3:
                    hdu.data = hdu.data[0:options.rows, :, :]
                    # hdu.add_checksum()
                if 'CHECKSUM' in hdu.header:
                    hdu.add_checksum()
            new_hdulist.writeto(os.path.join(options.output, b), output_verify='warn',
                                overwrite=options.overwrite, checksum=False)
    return 0


if __name__ == '__main__':
    sys.exit(main())
