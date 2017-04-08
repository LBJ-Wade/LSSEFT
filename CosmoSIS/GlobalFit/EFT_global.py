import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
from cosmosis.runtime.declare import declare_module
from astropy.io import ascii
from scipy import interpolate


class EFT_global(object):
    likes = section_names.likelihoods

    def __init__(self, my_config, my_name):
        self.mod_name = my_name

        z0 = self.import_z(my_config, my_name, 'fit_kmin_z0', 'fit_kmax_z0', 'data_z0_real', 'data_z0_multipole',
                           'theory_z0', 'real_sigma_z0', 'multipole_sigma_z0')
        z025 = self.import_z(my_config, my_name, 'fit_kmin_z025', 'fit_kmax_z025', 'data_z025_real',
                             'data_z025_multipole', 'theory_z025', 'real_sigma_z025', 'multipole_sigma_z025')
        z05 = self.import_z(my_config, my_name, 'fit_kmin_z05', 'fit_kmax_z05', 'data_z05_real',
                            'data_z05_multipole', 'theory_z05', 'real_sigma_z05', 'multipole_sigma_z05')
        z075 = self.import_z(my_config, my_name, 'fit_kmin_z075', 'fit_kmax_z075', 'data_z075_real',
                             'data_z075_multipole', 'theory_z075', 'real_sigma_z075', 'multipole_sigma_z075')
        z1 = self.import_z(my_config, my_name, 'fit_kmin_z1', 'fit_kmax_z1', 'data_z1_real', 'data_z1_multipole',
                           'theory_z1', 'real_sigma_z1', 'multipole_sigma_z1')

        self.payload = {'z0': z0, 'z025': z025, 'z05': z05, 'z075': z075, 'z1': z1}
        self.labels = ['z0', 'z025', 'z05', 'z075', 'z1']

    def import_z(self, my_config, my_name, fit_min, fit_max, data_real, data_mpole, theory, real_sigma, mpole_sigma):
        # get k-cuts for fit
        kmin = my_config[my_name, fit_min]
        kmax = my_config[my_name, fit_max]
        real_err = my_config[my_name, real_sigma]
        mpole_err = my_config[my_name, mpole_sigma]

        # import real-space and multipole data, and single file of theory data
        real_file = my_config[my_name, data_real]
        mpole_file = my_config[my_name, data_mpole]
        theory_file = my_config[my_name, theory]

        real = ascii.read(real_file)
        mpole = ascii.read(mpole_file)
        theory = ascii.read(theory_file)

        real = self.import_real(real, theory, kmin, kmax, real_err)
        P0 = self.import_mpole(mpole, theory, 'P0', kmin, kmax, mpole_err)
        P2 = self.import_mpole(mpole, theory, 'P2', kmin, kmax, mpole_err)
        P4 = self.import_mpole(mpole, theory, 'P4', kmin, kmax, mpole_err)

        dict = {'real': real, 'P0': P0, 'P2': P2, 'P4': P4}

        return dict

    def import_real(self, real, theory, kmin, kmax, err):
        # import theory k-values and predictions
        t_ks = theory['k']
        t_Pk = theory['dd']
        t_Z2_d = theory['dd_Z2_d']

        # import measured k-value and prediction
        d_ks = real['k/h']
        d_Pk = real['Pk']

        # cut data to fitted region
        lo_k = d_ks >= kmin
        hi_k = d_ks <= kmax
        mask = np.all([lo_k, hi_k], axis=0)

        d_ks_cut = d_ks[mask]
        d_Pk_cut = d_Pk[mask]

        # regrid theory predictions to cut data grid
        t_Pk_cut = self.regrid(t_ks, t_Pk, d_ks_cut)
        t_Z2_d_cut = self.regrid(t_ks, t_Z2_d, d_ks_cut)

        # precompute differences between SPT values and measured values
        Delta = t_Pk_cut - d_Pk_cut

        # estimate sigma at each k as a fixed fraction of the measured value
        Sigma = err * d_Pk_cut * err * d_Pk_cut

        dict = {'Delta': Delta, 'Sigma': Sigma, 'Z2_d': t_Z2_d_cut}

        return dict

    def import_mpole(self, real, theory, tag, kmin, kmax, err):
        # import theory k-values and predictions
        t_ks = theory['k']
        t_Pk = theory[tag]
        t_Z2_d = theory[tag + '_Z2_d']
        t_Z2_v = theory[tag + '_Z2_v']
        t_Z2_vd = theory[tag + '_Z2_vd']
        t_Z2_vv_A = theory[tag + '_Z2_vv_A']
        t_Z2_vv_B = theory[tag + '_Z2_vv_B']
        t_Z2_vvd = theory[tag + '_Z2_vvd']
        t_Z2_vvv = theory[tag + '_Z2_vvv']

        # import measured k-value and prediction
        d_ks = real['k/h']
        d_Pk = real[tag]

        # cut data to fitted region
        lo_k = d_ks >= kmin
        hi_k = d_ks <= kmax
        mask = np.all([lo_k, hi_k], axis=0)

        d_ks_cut = d_ks[mask]
        d_Pk_cut = d_Pk[mask]

        # regrid theory predictions to the cut data grid
        t_Pk_cut = self.regrid(t_ks, t_Pk, d_ks_cut)
        t_Z2_d_cut = self.regrid(t_ks, t_Z2_d, d_ks_cut)
        t_Z2_v_cut = self.regrid(t_ks, t_Z2_v, d_ks_cut)
        t_Z2_vd_cut = self.regrid(t_ks, t_Z2_vd, d_ks_cut)
        t_Z2_vv_A_cut = self.regrid(t_ks, t_Z2_vv_A, d_ks_cut)
        t_Z2_vv_B_cut = self.regrid(t_ks, t_Z2_vv_B, d_ks_cut)
        t_Z2_vvd_cut = self.regrid(t_ks, t_Z2_vvd, d_ks_cut)
        t_Z2_vvv_cut = self.regrid(t_ks, t_Z2_vvv, d_ks_cut)

        # precompute differences between SPT values and measured values
        Delta = t_Pk_cut - d_Pk_cut

        # estimate sigma at each k as a fixed fraction of the measured value
        Sigma = err * d_Pk_cut * err * d_Pk_cut

        dict = {'Delta': Delta, 'Sigma': Sigma, 'Z2_d': t_Z2_d_cut, 'Z2_v': t_Z2_v_cut, 'Z2_vd': t_Z2_vd_cut,
                'Z2_vv_A': t_Z2_vv_A_cut, 'Z2_vv_B': t_Z2_vv_B_cut, 'Z2_vvd': t_Z2_vvd_cut, 'Z2_vvv': t_Z2_vvv_cut}

        return dict

    def execute(self, block):
        Z2_d = block['EFT_counterterms', 'Z2_d']
        Z2_v = block['EFT_counterterms', 'Z2_v']
        Z2_vd = block['EFT_counterterms', 'Z2_vd']
        Z2_vv_A = block['EFT_counterterms', 'Z2_vv_A']
        Z2_vv_B = block['EFT_counterterms', 'Z2_vv_B']
        Z2_vvd = block['EFT_counterterms', 'Z2_vvd']
        Z2_vvv = block['EFT_counterterms', 'Z2_vvv']

        parameters = {'Z2_d': Z2_d, 'Z2_v': Z2_v, 'Z2_vd': Z2_vd, 'Z2_vv_A': Z2_vv_A, 'Z2_vv_B': Z2_vv_B,
                      'Z2_vvd': Z2_vvd, 'Z2_vvv': Z2_vvv}

        real_cterms = ['Z2_d']
        mpole_cterms = ['Z2_d', 'Z2_v', 'Z2_vd', 'Z2_vv_A', 'Z2_vv_B', 'Z2_vvd', 'Z2_vvv']

        L = 0

        for label in self.labels:
            z_bundle = self.payload[label]

            L += self.like(z_bundle['real'], parameters, real_cterms)
            L += self.like(z_bundle['P0'], parameters, mpole_cterms)
            L += self.like(z_bundle['P2'], parameters, mpole_cterms)
            L += self.like(z_bundle['P2'], parameters, mpole_cterms)

        block[self.likes, 'EFT_GLOBAL_LIKE'] = L

        return 0

    def like(self, dict, params, cterms):
        sum_list = [dict['Delta']]

        for label in cterms:
            Z = params[label]
            P_Z = dict[label]
            sum_list.append(Z * P_Z)

        Delta = np.sum(sum_list, axis=0)
        Sigma = dict['Sigma']

        return np.sum(- Delta * Delta / Sigma / 2.0)

    def cleanup(self):
        return 0

    def regrid(self, in_k, in_data, out_k):
        spline = interpolate.InterpolatedUnivariateSpline(in_k, in_data, ext='raise')
        return spline(out_k)


# register this module with the CosmoSIS core
declare_module(EFT_global)
