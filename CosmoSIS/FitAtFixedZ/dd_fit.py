import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
from cosmosis.runtime.declare import declare_module
from astropy.io import ascii
from scipy import interpolate


class dd_fit(object):

    likes = section_names.likelihoods

    def __init__(self, my_config, my_name):

        self.mod_name = my_name

        self.kmin = my_config[my_name, 'fit_kmin']
        self.kmax = my_config[my_name, 'fit_kmax']

        self.Delta_Pk, self.Z2d, self.sigma_Pk = self.import_realspace(my_config, my_name)

    def import_realspace(self, my_config, my_name):

        theory_file = my_config[my_name, 'realspace_theory']
        data_file = my_config[my_name, 'realspace_data']

        error = my_config[my_name, 'realspace_error']

        D_linear = my_config[my_name, 'D_linear']
        D_Zdelta = my_config[my_name, 'D_Zdelta']

        # compute rescaling factors for counterterms
        # (we normally want to report values for the c_i rather than the Z_i)
        Z2d_rescale = - D_Zdelta / (2 * D_linear * D_linear)

        theory = ascii.read(theory_file)
        data = ascii.read(data_file)

        # import predictions
        ks = theory['k']
        Pk = theory['dd']
        Z2d = theory['Z2_d']

        # import and mask data
        raw_data_ks = data['k/h']
        raw_data_Pk = data['Pk']

        lo_k_mask = (raw_data_ks >= self.kmin)
        hi_k_mask = (raw_data_ks <= self.kmax)
        mask = np.all([lo_k_mask, hi_k_mask], axis=0)

        data_ks = raw_data_ks[mask]
        data_Pk = raw_data_Pk[mask]

        # re-grid all theory quantities to the points at which we have data samples
        Pk = self.regrid(ks, Pk, data_ks)
        Z2d = self.regrid(ks, Z2d, data_ks) / Z2d_rescale

        # precompute difference between SPT value and measured value
        Delta_Pk = Pk - data_Pk

        # estimate sigma at each k as a fixed fraction of the power spectrum
        sigma_Pk = (error * data_Pk) * (error * data_Pk)

        return (Delta_Pk, Z2d, sigma_Pk)

    def execute(self, block):

        mu0 = block['dd_counterterms', 'c_mu0']

        Delta_Pk = self.Delta_Pk + mu0 * self.Z2d
        Dk = - Delta_Pk * Delta_Pk / self.sigma_Pk / 2.0

        like = np.sum(Dk)

        block[self.likes, 'DDFIT_LIKE'] = like

        return 0

    def cleanup(self):

        return 0

    def regrid(self, in_k, in_data, out_k):

        spline = interpolate.InterpolatedUnivariateSpline(in_k, in_data, ext='raise')
        return spline(out_k)


# register this module with the CosmoSIS core
declare_module(dd_fit)
