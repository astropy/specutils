import logging

import numpy as np
from astropy import units as u
from astropy.io import fits
from .spectrum1d import Spectrum1D

__all__ = ['XraySpectrum1D', 'AreaResponse', 'ResponseMatrix', 'arf_loader']

# For dealing with varied unit string choices
EV   = ['eV', 'ev']
KEV  = ['kev', 'keV']
ANGS = ['angs', 'Angs', 'Angstrom', 'angstrom', 'Angstroms', 'angstroms', 'A', 'a']

## An introduction to X-ray analysis can be found here:
## http://cxc.cfa.harvard.edu/xrayschool/talks/intro_xray_analysis.pdf

def _unit_parser(unit_string):
    """
    A function for matching unit strings to the correct Astropy unit
    """
    if unit_string in EV:
        return u.eV
    if unit_string in KEV:
        return u.keV
    if unit_string in ANGS:
        return u.angstrom

class XraySpectrum1D(Spectrum1D):
    """
    Spectrum properties specific to X-ray data

    Parameters
    ----------
    bin_lo : astropy.Quantity
        The left edges for bin values

    bin_hi : astropy.Quantity
        The right edges for bin values

    counts : astropy.Quantity
        Counts histogram for the X-ray spectrum

    exposure : astropy.Quantity
        Exposure time for the dataset

    arf : specutils.AreaResponse or string, default None
        Strings will be passed to AreaResponse.__init__
        All other input types are stored as arf attribute

    rmf : specutils.ResponseMatrix or string, default None
        Strings will be passed to ResponseMatrix.__init__
        All other input types are stored as rmf attribute

    rest_value : astropy.units.Quantity, default 0 Angstrom
        See Spectrum1D rest_value input

    Attributes
    ----------
    Inherits from Spectrum1D object.
    The following inputs are stored as additional attributes.

    bin_lo

    bin_hi

    counts

    exposure

    arf

    rmf
    """
    def __init__(self, bin_lo, bin_hi, counts, exposure,
                 arf=None, rmf=None, rest_value=0.0 * u.angstrom, **kwargs):

        bin_mid = 0.5 * (bin_lo + bin_hi)
        Spectrum1D.__init__(self, spectral_axis=bin_mid, flux=counts,
                            rest_value=rest_value, **kwargs)

        self.bin_lo = bin_lo
        self.bin_hi = bin_hi
        self.exposure = exposure
        self.assign_rmf(rmf)
        self.assign_arf(arf)

    # Convenience function for Xray people
    @property
    def counts(self):
        return self.flux

    def assign_arf(self, arf_inp, **kwargs):
        """
        Assign an area response file (ARF) object to the XraySpectrum1D object

        Input
        -----
        arf_inp : string
            File name for the area response file (FITS file)

        Returns
        -------
        Modifies the XraySpectrum1D.arf attribute
        """
        if isinstance(arf_inp, str):
            self.arf = arf_loader(arf_inp, **kwargs)
        else:
            self.arf = arf_inp

    def assign_rmf(self, rmf_inp):
        """
        Assign a response matrix file (RMF) object to the XraySpectrum1D object

        Input
        -----
        rmf_inp : string
            File name for the response matrix file (FITS file)

        Returns
        -------
        Modifies the XraySpectrum1D.rmf attribute
        """
        if isinstance(rmf_inp, str):
            self.rmf = ResponseMatrix(rmf_inp)
        else:
            self.rmf = rmf_inp
        return

    def apply_response(self, mflux, exposure=None):
        """
        Given a model flux spectrum, apply the response. In cases where the
        spectrum has both an area response file (ARF) and a response matrix
        file (RMF), apply both. Otherwise, apply whatever response is in RMF.

        The model flux spectrum *must* be created using the same units and
        bins as in the ARF (where the ARF exists)!

        Parameters
        ----------
        mflux : iterable
            A list or array with the model flux values in phot/s/cm^-2

        exposure : float, default None
            By default, the exposure stored in the ARF will be used to compute
            the total counts per bin over the effective observation time.
            In cases where this might be incorrect (e.g. for simulated spectra
            where the pha file might have a different exposure value than the
            ARF), this keyword provides the functionality to override the
            default behaviour and manually set the exposure time to use.

        store_model_counts : bool
            If True, the output will also be stored in self.model_counts

        Returns
        -------
        model_counts : numpy.ndarray
            The model spectrum in units of counts/bin

        If no ARF file exists, it will return the model flux after applying the RMF
        If no RMF file exists, it will return the model flux after applying the ARF (with a warning)
        If no ARF and no RMF, it will return the model flux spectrum (with a warning)
        """

        if self.arf is not None:
            mrate  = self.arf.apply_arf(mflux, exposure=exposure) # ct/bin with no RMF applied
        else:
            mrate = mflux # assumes ct/bin is the input with no RMF applied

        if self.rmf is not None:
            result = self.rmf.apply_rmf(mrate.value)
        else:
            result = mrate

        return result

## ----  Supporting response file objects

class ResponseMatrix(object):
    def __init__(self, filename=None, extension=None,
                 energ_lo=None, energ_hi=None, matrix=None, energ_unit=None,
                 offset=0.0, n_grp=np.array([]), f_chan=np.array([]),
                 n_chan=np.array([]), detchans=0):
        """
        Load a response matrix file (RMF) from a FITS file.

        filename : str
            The file name with the RMF file

        extension : str (default None)
            FITS file extension keyword, if the response matrix is stored
            under an extension other than "MATRIX" or "SPECRESP MATRIX"

        energ_lo : numpy.ndarray
            The lower edges of the energy bins

        energ_hi : numpy.ndarray
            The upper edges of the energy bins

        matrix : numpy.ndarray
            The redistribution matrix as a flattened 1D vector

        energ_unit : astropy.units.Unit
            Description of the energy units used

        offset : float

        n_grp : numpy.ndarray
            the Array with the number of channels in each
            channel set

        f_chan : numpy.ndarray
            The starting channel for each channel group;
            If an element i in n_grp > 1, then the resulting
            row entry in f_chan will be a list of length n_grp[i];
            otherwise it will be a single number

        n_chan : numpy.ndarray
            The number of channels in each channel group. The same
            logic as for f_chan applies

        detchans : int
            The number of channels in the detector
        """
        self.filename = filename
        self.offset = offset
        self.n_grp = n_grp
        self.f_chan = f_chan
        self.n_chan = n_chan
        self.matrix = matrix
        self.energ_lo = energ_lo
        self.energ_hi = energ_hi
        self.energ_unit = energ_unit
        self.detchans = detchans
        if filename is not None:
            self._load_rmf(filename, extension=extension)

    def _load_rmf(self, filename, extension=None):
        # open the FITS file and extract the MATRIX extension
        # which contains the redistribution matrix and
        # anxillary information
        hdulist = fits.open(filename)

        # get all the extension names
        extnames = np.array([h.name for h in hdulist])

        # figure out the right extension to use
        if extension is not None:
            h = hdulist[extension]

        elif "MATRIX" in extnames:
            h = hdulist["MATRIX"]

        elif "SPECRESP MATRIX" in extnames:
            h = hdulist["SPECRESP MATRIX"]

        else:
            print("Cannot find common FITS file extension for response matrix values")
            print("Please set the `extension` keyword")
            return

        data = h.data
        hdr = h.header
        hdulist.close()

        # extract + store the attributes described in the docstring
        n_grp = np.array(data.field("N_GRP"))
        f_chan = np.array(data.field('F_CHAN'))
        n_chan = np.array(data.field("N_CHAN"))
        matrix = np.array(data.field("MATRIX"))

        self.energ_lo = np.array(data.field("ENERG_LO"))
        self.energ_hi = np.array(data.field("ENERG_HI"))
        self.energ_unit = u.Unit(data.columns["ENERG_LO"].unit)
        self.detchans = hdr["DETCHANS"]
        self.offset = self.__get_tlmin(h)

        # flatten the variable-length arrays
        self.n_grp, self.f_chan, self.n_chan, self.matrix = \
                self._flatten_arrays(n_grp, f_chan, n_chan, matrix)

        return

    def __get_tlmin(self, h):
        """
        Get the tlmin keyword for `F_CHAN`.

        Parameters
        ----------
        h : an astropy.io.fits.hdu.table.BinTableHDU object
            The extension containing the `F_CHAN` column

        Returns
        -------
        tlmin : int
            The tlmin keyword
        """
        # get the header
        hdr = h.header
        # get the keys of all
        keys = np.array(list(hdr.keys()))

        # find the place where the tlmin keyword is defined
        t = np.array(["TLMIN" in k for k in keys])

        # get the index of the TLMIN keyword
        tlmin_idx = np.hstack(np.where(t))[0]

        # get the corresponding value
        tlmin = np.int(list(hdr.items())[tlmin_idx][1])

        return tlmin

    def _flatten_arrays(self, n_grp, f_chan, n_chan, matrix):

        if not len(n_grp) == len(f_chan) == len(n_chan) == len(matrix):
            raise ValueError("Arrays must be of same length!")

        # find all non-zero groups
        nz_idx = (n_grp > 0)

        # stack all non-zero rows in the matrix
        matrix_flat = np.hstack(matrix[nz_idx])

        # stack all nonzero rows in n_chan and f_chan
        #n_chan_flat = np.hstack(n_chan[nz_idx])
        #f_chan_flat = np.hstack(f_chan[nz_idx])

        # some matrices actually have more elements
        # than groups in `n_grp`, so we'll only pick out
        # those values that have a correspondence in
        # n_grp
        f_chan_new = []
        n_chan_new = []
        for i,t in enumerate(nz_idx):
            if t:
                n = n_grp[i]
                f = f_chan[i]
                nc = n_chan[i]
                if np.size(f) == 1:
                    f_chan_new.append(f)
                    n_chan_new.append(nc)
                else:
                    f_chan_new.append(f[:n])
                    n_chan_new.append(nc[:n])

        n_chan_flat = np.hstack(n_chan_new)
        f_chan_flat = np.hstack(f_chan_new)

        # if n_chan is zero, we'll remove those as well.
        nz_idx2 = (n_chan_flat > 0)
        n_chan_flat = n_chan_flat[nz_idx2]
        f_chan_flat = f_chan_flat[nz_idx2]

        return n_grp, f_chan_flat, n_chan_flat, matrix_flat

    def apply_rmf(self, model):
        """
        Fold the spectrum through the redistribution matrix.

        The redistribution matrix is saved as a flattened 1-dimensional
        vector to save space. In reality, for each entry in the flux
        vector, there exists one or more sets of channels that this
        flux is redistributed into. The additional arrays n_grp,
        f_chan and n_chan store this information:

            * n_group stores the number of channel groups for each
              energy bin

            * f_chan stores the *first channel* that each channel
              for each channel set

            * n_chan stores the number of channels in each channel
              set

        As a result, for a given energy bin i, we need to look up the
        number of channel sets in n_grp for that energy bin. We
        then need to loop over the number of channel sets. For each
        channel set, we look up the first channel into which flux
        will be distributed as well as the number of channels in the
        group. We then need to also loop over the these channels and
        actually use the corresponding elements in the redistribution
        matrix to redistribute the photon flux into channels.

        All of this is basically a big bookkeeping exercise in making
        sure to get the indices right.

        Parameters
        ----------
        model : numpy.ndarray
            The (model) spectrum to be folded

        Returns
        -------
        counts : numpy.ndarray
            The (model) spectrum after folding, in
            counts/s/channel

        """
        # get the number of channels in the data
        nchannels = model.shape[0]

        # an empty array for the output counts
        counts = np.zeros(nchannels)

        # index for n_chan and f_chan incrementation
        k = 0

        # index for the response matrix incrementation
        resp_idx = 0

        # loop over all channels
        for i in range(nchannels):

            # this is the current bin in the flux spectrum to
            # be folded
            source_bin_i = model[i]

            # get the current number of groups
            current_num_groups = self.n_grp[i]

            # loop over the current number of groups
            for j in range(current_num_groups):

                current_num_chans = int(self.n_chan[k])

                if current_num_chans == 0:
                    k += 1
                    resp_idx += current_num_chans
                    continue


                else:
                    # get the right index for the start of the counts array
                    # to put the data into
                    counts_idx = int(self.f_chan[k] - self.offset)
                    # this is the current number of channels to use

                    k += 1
                    # add the flux to the subarray of the counts array that starts with
                    # counts_idx and runs over current_num_chans channels
                    counts[counts_idx:counts_idx +
                                      current_num_chans] += self.matrix[resp_idx:resp_idx +
                                                                                 current_num_chans] * \
                                                                np.float(source_bin_i)
                    # iterate the response index for next round
                    resp_idx += current_num_chans


        return counts[:self.detchans]


class AreaResponse(object):
    def __init__(self, e_low, e_high, eff_area,
                 exposure=None, fracexpo=1.0, filename=None):
        """
        Load an area response file (ARF) from a FITS file.

        Attributes
        ----------
        e_low : astropy.units.Quantity
            The lower edges of the energy bins

        e_high : astropy.units.Quantity
            The upper edges of the energy bins

        eff_area : astropy.units.Quantity
            Description of the energy dependent telescope effective area

        exposure :
            Average exposure time for the dataset
            (which takes telescope dithering into account).

        fracexpo : float or numpy.ndarray
            Fractional exposure time for the spectrum
            (sometimes constant, sometimes dependent on spectral channel).
            These values are stored for reference; generally, they are already
            accounted for in the eff_area array.

        filename : str
            Stores the name of the file from which the ARF was loaded.

        e_mid : Property that returns the middle of each energy bin
        """
        self.filename = filename
        self.e_low    = e_low
        self.e_high   = e_high
        self.eff_area = eff_area
        self.fracexpo = fracexpo
        self.exposure = exposure

    @property
    def e_mid(self):
        return 0.5 * (self.e_low + self.e_high)

    def apply_arf(self, model, exposure=None):
        """
        Fold the spectrum through the area response file (ARF).
        The ARF is a single vector encoding the effective area information
        about the detector. A such, applying the ARF is a simple
        multiplication with the input spectrum.

        Parameters
        ----------
        model : numpy.ndarray
            The model spectrum to which the arf will be applied

        exposure : float, default None
            Value for the exposure time. By default, `apply_arf` will use the
            exposure keyword from the ARF file. If this exposure time is not
            correct (for example when simulated spectra use a different exposure
            time and the ARF from a real observation), one can override the
            default exposure by setting the `exposure` keyword to the correct
            value.

        Returns
        -------
        s_arf : numpy.ndarray
            The (model) spectrum after folding, in
            counts/s/channel
        """
        assert model.shape[0] == self.eff_area.shape[0], "The input spectrum must " \
                                                      "be of same size as the " \
                                                      "ARF array."
        if exposure is None:
            return model * self.eff_area * self.exposure
        else:
            return model * self.eff_area * exposure

## ----- Convenience functions for loading ARFs and RMFs
def arf_loader(filename, block='SPECRESP'):
    """
    Load a AreaResponse object from FITS file.

    Parameters
    ----------
    filename : str
        Path to the FITS file

    block : str
        FITS file block keyword, if the spectral response is stored
        under an extension other than "SPECRESP"

    Returns
    -------
    AreaResponse
        The ARF that is represented by the FITS file
    """
    # open the FITS file and extract the SPECRESP block
    # which contains the spectral response for the telescope
    hdulist = fits.open(filename)

    # Get the data from the appropriate extension
    blocknames = np.array([h.name for h in hdulist])
    assert block in blocknames, "Could not fine {} block in FITS file".format(block)
    h = hdulist[block]

    data = h.data
    hdr = h.header
    hdulist.close()

    e_unit = u.Unit(data.columns["ENERG_LO"].unit)
    e_low  = np.array(data.field("ENERG_LO")) * e_unit
    e_high = np.array(data.field("ENERG_HI")) * e_unit

    area_unit = u.Unit(data.columns['SPECRESP'].unit) * u.ct # usually ct cm^2 / phot
    specresp  = np.array(data.field("SPECRESP")) * area_unit

    exposure = None
    if "EXPOSURE" in list(hdr.keys()):
        exposure = hdr["EXPOSURE"]

    fracexpo = 1.0
    if "FRACEXPO" in data.columns.names:
        fracexpo = data["FRACEXPO"]

    return AreaResponse(e_low, e_high, specresp,
                        exposure=exposure, fracexpo=fracexpo, filename=filename)
