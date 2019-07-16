import math
from astropy import units as u
from specutils import Spectrum1D


class TempSpec:
    """
    Temporary class to assemble Spectrum1D parameters from ascii files
    """
    def __init__(self, filename):
        self.file_data_wave = []
        self.file_data_flux = []
        self.file_data_weight = []

        self.filename = filename
        self.load_file(filename)

    def load_file(self, filename):
        first_pass = True

        for line in open(filename, 'r'):
            if first_pass:
                first_pass = False
                continue
            item = line.rstrip()
            floats = list(map(float, item.split()))
            self.add_floats(floats)

    def add_floats(self, floats):
        self.file_data_wave.append(floats[0])
        self.file_data_flux.append(floats[1])
        if len(floats) > 2:
            self.file_data_weight.append(floats[2])

class TemplateMatching:
    """
    Takes in two files (first with observed wavelengths, second with template/model wavelengths) and finds the
    chi square of the wavelengths the two spectra have in commmon
    """

    def __init__(self, observed_filename, template_filename):
        self.weightj = 0.0
        self.chi_square = 0.0

        # load files into temporary objects to be assembled to later be converted into Spectrum1D objects
        self.observed = TempSpec(observed_filename)
        self.template = TempSpec(template_filename)

        lamb_observed = self.observed.file_data_wave * u.AA  # doctest: +REMOTE_DATA
        flux_observed = self.observed.file_data_flux * u.Unit('erg cm-2 s-1 AA-1')  # doctest: +REMOTE_DATA
        self.spec_observed = Spectrum1D(spectral_axis=lamb_observed, flux=flux_observed)  # doctest: +REMOTE_DATA

        lamb_template = self.template.file_data_wave * u.AA  # doctest: +REMOTE_DATA
        flux_template = self.template.file_data_flux * u.Unit('erg cm-2 s-1 AA-1')  # doctest: +REMOTE_DATA
        self.spec_template = Spectrum1D(spectral_axis=lamb_template, flux=flux_template)  # doctest: +REMOTE_DATA

    def get_weight_numerator(self, legac_flux, model_flux, weight):
        """
        Assemble weightj numerator during the first pass through the spectra
        """
        sigma = self.get_sigma_from_weight(weight)

        # Use flux.value since there a problems with squaring certain units
        numerator = (legac_flux.value * model_flux.value) / (math.pow(sigma, 2))
        return numerator

    def get_weight_denominator(self, model_flux, weight):
        """
        Assemble weightj numerator during the first pass through the spectra
        """
        sigma = self.get_sigma_from_weight(weight)

        # Use flux.value since there a problems with squaring certain units
        denominator = math.pow((model_flux.value / sigma), 2)
        return denominator

    def get_chi_square(self, legac_flux, weightj, model_flux, weight):
        sigma = self.get_sigma_from_weight(weight)

        # Use flux.value since there a problems with squaring certain units
        chi_square = math.pow(((legac_flux.value - (weightj * model_flux.value)) / sigma), 2)
        return chi_square

    def get_sigma_from_weight(self, weight):
        """
        Convert from weight to sigma
        """

        if not weight == 0.0:
            return 1 / math.sqrt(weight)
        return 0.0

    def run_chi_square(self):
        """
        Runs the chi square formula on an observed spectrum and a template spectrum and returns the result
        """
        observed_wavelength = self.spec_observed.wavelength
        template_wavelength = self.spec_template.wavelength

        model_min = 0
        model_end = len(template_wavelength)
        observed_end = len(observed_wavelength)

        wj_numerator = 0.0
        wj_denominator = 0.0

        for wave_o in range(0, observed_end):
            if self.spec_observed.file_data_weight[wave_o] == 0.0:
                continue
            for wave_t in range(model_min, model_end):
                if observed_wavelength[wave_o] > template_wavelength[wave_t] - (0.1 * template_wavelength.unit) and \
                        observed_wavelength[wave_o] < template_wavelength[wave_t] + (0.1 * template_wavelength.unit):

                    wj_numerator += self.get_weight_numerator(self.spec_observed.flux[wave_o],
                                                              self.spec_template.flux[wave_t],
                                                              self.observed.file_data_weight[wave_o])
                    wj_denominator += self.get_weight_denominator(self.spec_template.flux[wave_t],
                                                                  self.observed.file_data_weight[wave_o])
                    model_min = wave_t
                    break
                elif observed_wavelength[wave_o] < template_wavelength[wave_t] + (0.1 * template_wavelength.unit):
                    model_min = wave_t
                    break

        self.weightj = wj_numerator / wj_denominator
        model_min = 0

        for wave_o in range(0, observed_end):
            if self.observed.file_data_weight[wave_o] == 0.0:
                continue
            for wave_t in range(model_min, model_end):
                if observed_wavelength[wave_o] > template_wavelength[wave_t] - (0.1 * template_wavelength.unit) and \
                        observed_wavelength[wave_o] < template_wavelength[wave_t] + (0.1 * template_wavelength.unit):
                    #             print(observed_wavelength[b])

                    self.chi_square += self.get_chi_square(observed_wavelength[wave_o], self.weightj,
                                                      template_wavelength[wave_t],
                                                      self.observed.file_data_weight[wave_o])
                    model_min = wave_t
                    break
                elif observed_wavelength[wave_o] < template_wavelength[wave_t] + (0.1 * template_wavelength.unit):
                    model_min = wave_t
                    break
        self.chi_square = self.chi_square * (-0.5)
        return self.chi_square
