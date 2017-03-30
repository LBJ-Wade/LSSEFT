import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
from cosmosis.runtime.declare import declare_module
from astropy.io import ascii
from scipy import interpolate


class rsd_fit(object):

    likes = section_names.likelihoods

    def __init__(self, my_config, my_name):

        self.mod_name = my_name

        self.kmin = my_config[my_name, 'fit_kmin']
        self.kmax = my_config[my_name, 'fit_kmax']

        # need growth factor to compute mu6 coefficient
        self.f = my_config[my_name, 'f']

        # import theory predictions and measured values for multipoles
        P0_group, P2_group, P4_group = self.import_multipoles(my_config, my_name)

        self.Delta_P0, self.P0_mu0, self.P0_mu2, self.P0_mu4, self.P0_mu6, self.sigma_P0 = P0_group
        self.Delta_P2, self.P2_mu0, self.P2_mu2, self.P2_mu4, self.P2_mu6, self.sigma_P2 = P2_group
        self.Delta_P4, self.P4_mu0, self.P4_mu2, self.P4_mu4, self.P4_mu6, self.sigma_P4 = P4_group

        # determine whether we are using the realspace power spectrum to
        # anchor the fit
        self.use_realspace = my_config[my_name, 'use_realspace']

        if self.use_realspace:
            print 'Using realspace power spectrum to assist fit'
            self.Delta_Pk, self.Z2d, self.sigma_Pk = self.import_realspace(my_config, my_name)

    def import_multipoles(self, my_config, my_name):

        theory_file = my_config[my_name, 'multipole_theory']
        data_file = my_config[my_name, 'multipole_data']

        D_linear = my_config[my_name, 'D_linear']
        D_mu0 = my_config[my_name, 'D_mu0']
        D_mu2 = my_config[my_name, 'D_mu2']
        D_mu4 = my_config[my_name, 'D_mu4']
        D_mu6 = my_config[my_name, 'D_mu6']

        # compute rescaling factors for counterterms
        # (we normally want to report values for the c_i rather than the Z_i)
        mu0_rescale = - D_mu0 / (2 * D_linear * D_linear)
        mu2_rescale = - D_mu2 / (2 * D_linear * D_linear)
        mu4_rescale = - D_mu4 / (2 * D_linear * D_linear)
        mu6_rescale = - D_mu6 / (2 * D_linear * D_linear)

        theory = ascii.read(theory_file)
        data = ascii.read(data_file)

        ks = theory['k']
        P0 = theory['P0']
        P2 = theory['P2']
        P4 = theory['P4']

        P0_mu0 = theory['P0_mu0']
        P0_mu2 = theory['P0_mu2']
        P0_mu4 = theory['P0_mu4']
        P0_mu6 = theory['P0_mu6']

        P2_mu0 = theory['P2_mu0']
        P2_mu2 = theory['P2_mu2']
        P2_mu4 = theory['P2_mu4']
        P2_mu6 = theory['P2_mu6']

        P4_mu0 = theory['P4_mu0']
        P4_mu2 = theory['P4_mu2']
        P4_mu4 = theory['P4_mu4']
        P4_mu6 = theory['P4_mu6']

        # import and mask data
        raw_data_ks = data['k/h']
        raw_data_P0 = data['P0']
        raw_data_P2 = data['P2']
        raw_data_P4 = data['P4']

        lo_k_mask = (raw_data_ks >= self.kmin)
        hi_k_mask = (raw_data_ks <= self.kmax)
        mask = np.all([lo_k_mask, hi_k_mask], axis=0)

        data_ks = raw_data_ks[mask]
        data_P0 = raw_data_P0[mask]
        data_P2 = raw_data_P2[mask]
        data_P4 = raw_data_P4[mask]

        # re-grid all theory quantities to the points at which we have data samples
        P0 = self.regrid(ks, P0, data_ks)
        P2 = self.regrid(ks, P2, data_ks)
        P4 = self.regrid(ks, P4, data_ks)

        P0_mu0 = self.regrid(ks, P0_mu0, data_ks) / mu0_rescale
        P0_mu2 = self.regrid(ks, P0_mu2, data_ks) / mu2_rescale
        P0_mu4 = self.regrid(ks, P0_mu4, data_ks) / mu4_rescale
        P0_mu6 = self.regrid(ks, P0_mu6, data_ks) / mu6_rescale

        P2_mu0 = self.regrid(ks, P2_mu0, data_ks) / mu0_rescale
        P2_mu2 = self.regrid(ks, P2_mu2, data_ks) / mu2_rescale
        P2_mu4 = self.regrid(ks, P2_mu4, data_ks) / mu4_rescale
        P2_mu6 = self.regrid(ks, P2_mu6, data_ks) / mu6_rescale

        P4_mu0 = self.regrid(ks, P4_mu0, data_ks) / mu0_rescale
        P4_mu2 = self.regrid(ks, P4_mu2, data_ks) / mu2_rescale
        P4_mu4 = self.regrid(ks, P4_mu4, data_ks) / mu4_rescale
        P4_mu6 = self.regrid(ks, P4_mu6, data_ks) / mu6_rescale

        # precompute differences between SPT values and measured values
        Delta_P0 = P0 - data_P0
        Delta_P2 = P2 - data_P2
        Delta_P4 = P4 - data_P4

        # estimate sigma at each k as 20% of the power spectrum
        sigma_P0 = (0.20 * data_P0) * (0.20 * data_P0)
        sigma_P2 = (0.20 * data_P2) * (0.20 * data_P2)
        sigma_P4 = (0.20 * data_P4) * (0.20 * data_P4)

        P0_group = (Delta_P0, P0_mu0, P0_mu2, P0_mu4, P0_mu6, sigma_P0)
        P2_group = (Delta_P2, P2_mu0, P2_mu2, P2_mu4, P2_mu6, sigma_P2)
        P4_group = (Delta_P4, P4_mu0, P4_mu2, P4_mu4, P4_mu6, sigma_P4)

        return (P0_group, P2_group, P4_group)

    def import_realspace(self, my_config, my_name):

        theory_file = my_config[my_name, 'realspace_theory']
        data_file = my_config[my_name, 'realspace_data']

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

        # estimate sigma at each k as 5% of the power spectrum
        sigma_Pk = (0.05 * data_Pk) * (0.05 * data_Pk)

        return (Delta_Pk, Z2d, sigma_Pk)

    def execute(self, block):

        mu0 = block['rsd_counterterms', 'c_mu0']
        mu2 = block['rsd_counterterms', 'c_mu2']
        mu4 = block['rsd_counterterms', 'c_mu4']
        mu6 = self.f * self.f * self.f * mu0 - self.f * self.f * mu2 + self.f * mu4

        Delta_P0 = self.Delta_P0 + mu0 * self.P0_mu0 + mu2 * self.P0_mu2 + mu4 * self.P0_mu4 + mu6 * self.P0_mu6
        Delta_P2 = self.Delta_P2 + mu0 * self.P2_mu0 + mu2 * self.P2_mu2 + mu4 * self.P2_mu4 + mu6 * self.P2_mu6
        Delta_P4 = self.Delta_P4 + mu0 * self.P4_mu0 + mu2 * self.P4_mu2 + mu4 * self.P4_mu4 + mu6 * self.P4_mu6

        D0 = - Delta_P0 * Delta_P0 / self.sigma_P0 / 2.0
        D2 = - Delta_P2 * Delta_P2 / self.sigma_P2 / 2.0
        D4 = - Delta_P4 * Delta_P4 / self.sigma_P4 / 2.0

        like = np.sum(D0) + np.sum(D2) + np.sum(D4)

        if self.use_realspace:
            Delta_Pk = self.Delta_Pk + mu0 * self.Z2d
            Dk = - Delta_Pk * Delta_Pk / self.sigma_Pk / 2.0

            like += np.sum(Dk)

        block[self.likes, 'RSDFIT_LIKE'] = like

        return 0

    def cleanup(self):

        return 0

    def regrid(self, in_k, in_data, out_k):

        spline = interpolate.InterpolatedUnivariateSpline(in_k, in_data, ext='raise')
        return spline(out_k)


# register this module with the CosmoSIS core
declare_module(rsd_fit)
