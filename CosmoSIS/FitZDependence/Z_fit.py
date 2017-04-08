import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
from cosmosis.runtime.declare import declare_module
from astropy.io import ascii
from scipy import interpolate


class Z_fit(object):

    likes = section_names.likelihoods

    def __init__(self, my_config, my_name):

        self.mod_name = my_name

        data_file = my_config[my_name, 'data']
        mu0_file = my_config[my_name, 'mu0']
        mu2_file = my_config[my_name, 'mu2']
        mu4_file = my_config[my_name, 'mu4']

        data = ascii.read(data_file)
        mu0 = ascii.read(mu0_file)
        mu2 = ascii.read(mu2_file)
        mu4 = ascii.read(mu4_file)

        self.data_c0 = data['c_0']
        self.data_c2 = data['c_2']
        self.data_c4 = data['c_4']

        # estimate error in each measured value as 5%
        self.sigma_c0 = 0.05*self.data_c0 * 0.05*self.data_c0
        self.sigma_c2 = 0.05*self.data_c2 * 0.05*self.data_c2
        self.sigma_c4 = 0.05*self.data_c4 * 0.05*self.data_c4

        self.mu0_Z2_d = mu0['Z2_d']
        self.mu0_Z2_v = mu0['Z2_v']
        self.mu0_Z2_vd = mu0['Z2_vd']
        self.mu0_Z2_vv_A = mu0['Z2_vv_A']
        self.mu0_Z2_vv_B = mu0['Z2_vv_B']
        self.mu0_Z2_vvd = mu0['Z2_vvd']
        self.mu0_Z2_vvv = mu0['Z2_vvv']

        self.mu2_Z2_d = mu2['Z2_d']
        self.mu2_Z2_v = mu2['Z2_v']
        self.mu2_Z2_vd = mu2['Z2_vd']
        self.mu2_Z2_vv_A = mu2['Z2_vv_A']
        self.mu2_Z2_vv_B = mu2['Z2_vv_B']
        self.mu2_Z2_vvd = mu2['Z2_vvd']
        self.mu2_Z2_vvv = mu2['Z2_vvv']

        self.mu4_Z2_d = mu4['Z2_d']
        self.mu4_Z2_v = mu4['Z2_v']
        self.mu4_Z2_vd = mu4['Z2_vd']
        self.mu4_Z2_vv_A = mu4['Z2_vv_A']
        self.mu4_Z2_vv_B = mu4['Z2_vv_B']
        self.mu4_Z2_vvd = mu4['Z2_vvd']
        self.mu4_Z2_vvv = mu4['Z2_vvv']

    def execute(self, block):

        Z2_d = block['zdep_counterterms', 'Z2_d']
        Z2_v = block['zdep_counterterms', 'Z2_v']
        Z2_vd = block['zdep_counterterms', 'Z2_vd']
        Z2_vv_A = block['zdep_counterterms', 'Z2_vv_A']
        Z2_vv_B = block['zdep_counterterms', 'Z2_vv_B']
        Z2_vvd = block['zdep_counterterms', 'Z2_vvd']
        Z2_vvv = block['zdep_counterterms', 'Z2_vvv']

        c0 = Z2_d * self.mu0_Z2_d + Z2_v * self.mu0_Z2_v + Z2_vd * self.mu0_Z2_vd + Z2_vv_A * self.mu0_Z2_vv_A + Z2_vv_B * self.mu0_Z2_vv_B + Z2_vvd * self.mu0_Z2_vvd + Z2_vvv * self.mu0_Z2_vvv
        c2 = Z2_d * self.mu2_Z2_d + Z2_v * self.mu2_Z2_v + Z2_vd * self.mu2_Z2_vd + Z2_vv_A * self.mu2_Z2_vv_A + Z2_vv_B * self.mu2_Z2_vv_B + Z2_vvd * self.mu2_Z2_vvd + Z2_vvv * self.mu2_Z2_vvv
        c4 = Z2_d * self.mu4_Z2_d + Z2_v * self.mu4_Z2_v + Z2_vd * self.mu4_Z2_vd + Z2_vv_A * self.mu4_Z2_vv_A + Z2_vv_B * self.mu4_Z2_vv_B + Z2_vvd * self.mu4_Z2_vvd + Z2_vvv * self.mu4_Z2_vvv

        Delta_c0 = c0 - self.data_c0
        Delta_c2 = c2 - self.data_c2
        Delta_c4 = c4 - self.data_c4

        lik_c0 = - Delta_c0 * Delta_c0 / self.sigma_c0 / 2.0
        lik_c2 = - Delta_c2 * Delta_c2 / self.sigma_c2 / 2.0
        lik_c4 = - Delta_c4 * Delta_c4 / self.sigma_c4 / 2.0

        lik = np.sum(lik_c0) + np.sum(lik_c2) + np.sum(lik_c4)

        block[self.likes, 'Z_FIT_LIKE'] = lik

    def cleanup(self):

        return 0


# register this module with the CosmoSIS core
declare_module(Z_fit)
