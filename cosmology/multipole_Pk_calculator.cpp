//
// Created by David Seery on 18/11/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#include <cmath>

#include "multipole_Pk_calculator.h"

#include "cuba.h"

namespace multipole_Pk_calculator_impl
  {

    constexpr unsigned int dimensions            = 2;   // no point doing integrals over phi, because the integrands don't depend on it
    constexpr unsigned int components            = 1;
    constexpr unsigned int points_per_invocation = 1;
    
    constexpr unsigned int verbosity_none        = 0;
    constexpr unsigned int verbosity_reasonable  = 1;
    constexpr unsigned int verbosity_progress    = 2;
    constexpr unsigned int verbosity_subregions  = 3;
    
    constexpr unsigned int samples_all           = 0;
    constexpr unsigned int samples_last          = 4;
    
    constexpr unsigned int min_eval              = 0;
    constexpr unsigned int max_eval              = 20000000;
    
    constexpr unsigned int ngiven                = 0;
    constexpr unsigned int ldxgiven              = 0;
    constexpr unsigned int nextra                = 0;
    
    constexpr unsigned int pcores                = 10000;   // matches default Cuba value
    
    constexpr unsigned int cuhre_key             = 13;      // degree-13 only available in 2-dimensions
    constexpr unsigned int cuhre_key_1d          = 9;
    constexpr unsigned int divonne_key1          = 47;
    constexpr unsigned int divonne_key2          = 13;      // degree-13 only available in 2-dimensions
    constexpr unsigned int divonne_key3          = 1;
    constexpr unsigned int divonne_maxpass       = 5;
    constexpr unsigned int divonne_border        = 0;
    constexpr double       divonne_maxchisq      = 10.0;
    constexpr double       divonne_minchisq      = 0.25;
    
    
    class integrand_data
      {
      
      public:
        
        integrand_data(const Mpc_units::energy& _k, const Mpc_units::energy& UV, const Mpc_units::energy& IR,
                       const tree_power_spectrum& _Pk)
          : k(_k),
            UV_cutoff(UV),
            IR_cutoff(IR),
            Pk(_Pk),
            jacobian(UV_cutoff - IR_cutoff),
            q_range(UV_cutoff - IR_cutoff),
            k_sq(k*k)
          {
          }
        
        const Mpc_units::energy& k;
        const Mpc_units::energy& UV_cutoff;
        const Mpc_units::energy& IR_cutoff;
        const tree_power_spectrum& Pk;
        
        Mpc_units::energy  jacobian;
        Mpc_units::energy  q_range;
        Mpc_units::energy2 k_sq;
      };
    
    
    struct mu_to_ell0
      {
        
        double operator()(mu_power n)
          {
            switch(n)
              {
                case mu_power::mu0:
                  return 2.0;
                case mu_power::mu2:
                  return 2.0 / 3.0;
                case mu_power::mu4:
                  return 2.0 / 5.0;
                case mu_power::mu6:
                  return 2.0 / 7.0;
                case mu_power::mu8:
                  return 2.0 / 9.0;
              }
          }
        
      };
    
    struct mu_to_ell0_rs
      {
        
        mu_to_ell0_rs(double _A, double _B)
          : A(_A),
            B(_B)
          {
          }
        
        double operator()(mu_power n)
          {
            switch(n)
              {
                case mu_power::mu0:
                  return (std::sqrt(M_PI) * std::erf(std::sqrt(A))) / (std::sqrt(A) * std::exp(B));
                case mu_power::mu2:
                  return
                    (std::exp(-A - B) * (-2. * std::sqrt(A) + std::exp(A) * std::sqrt(M_PI) * std::erf(std::sqrt(A)))) /
                    (2. * std::pow(A, 1.5));
                case mu_power::mu4:
                  return (std::exp(-A - B) * (-2. * std::sqrt(A) * (3. + 2. * A) +
                                              3. * std::exp(A) * std::sqrt(M_PI) * std::erf(std::sqrt(A)))) /
                         (4. * std::pow(A, 2.5));
                case mu_power::mu6:
                  return (std::exp(-A - B) * (-2. * std::sqrt(A) * (15. + 2. * A * (5. + 2. * A)) +
                                              15. * std::exp(A) * std::sqrt(M_PI) * std::erf(std::sqrt(A)))) /
                         (8. * std::pow(A, 3.5));
                case mu_power::mu8:
                  return (std::exp(-A - B) * (-2. * std::sqrt(A) * (105. + 2. * A * (35. + 2. * A * (7. + 2. * A))) +
                                              105. * std::exp(A) * std::sqrt(M_PI) * std::erf(std::sqrt(A)))) /
                         (16. * std::pow(A, 4.5));
              }
          }
        
        double A;
        double B;
        
      };

    struct mu_to_ell2
      {
        
        double operator()(mu_power n)
          {
            switch(n)
              {
                case mu_power::mu0:
                  return 0.0;
                case mu_power::mu2:
                  return 4.0 / 15.0;
                case mu_power::mu4:
                  return 8.0 / 35.0;
                case mu_power::mu6:
                  return 4.0 / 21.0;
                case mu_power::mu8:
                  return 16.0 / 99.0;
              }
          }
        
      };
    
    struct mu_to_ell2_rs
      {
        
        mu_to_ell2_rs(double _A, double _B)
          : A(_A),
            B(_B)
          {
          }
        
        double operator()(mu_power n)
          {
            switch(n)
              {
                case mu_power::mu0:
                  return (std::exp(-A - B) *
                          (-6. * std::sqrt(A) - (-3. + 2. * A) * std::exp(A) * std::sqrt(M_PI) * std::erf(std::sqrt(A)))) /
                         (4. * std::pow(A, 1.5));
                case mu_power::mu2:
                  return (std::exp(-A - B) * (-2. * std::sqrt(A) * (9. + 4. * A) -
                                              (-9. + 2. * A) * std::exp(A) * std::sqrt(M_PI) * std::erf(std::sqrt(A)))) /
                         (8. * std::pow(A, 2.5));
                case mu_power::mu4:
                  return (std::exp(-A - B) * (-2. * std::sqrt(A) * (45. + 8. * A * (3. + A)) -
                                              3. * (-15. + 2. * A) * std::exp(A) * std::sqrt(M_PI) *
                                              std::erf(std::sqrt(A)))) / (16. * std::pow(A, 3.5));
                case mu_power::mu6:
                  return (std::exp(-A - B) * (-2. * std::sqrt(A) * (315. + 4. * A * (45. + 4. * A * (4. + A))) -
                                              15. * (-21. + 2. * A) * std::exp(A) * std::sqrt(M_PI) *
                                              std::erf(std::sqrt(A)))) / (32. * std::pow(A, 4.5));
                case mu_power::mu8:
                  return (std::exp(-A - B) * (-2. * std::sqrt(A) * (2835. + 8. * A * (210. + A * (77. + 4. * A * (5. + A)))) -
                                              105. * (-27. + 2. * A) * std::exp(A) * std::sqrt(M_PI) *
                                              std::erf(std::sqrt(A)))) / (64. * std::pow(A, 5.5));
              }
          }
        
        double A;
        double B;
        
      };
    
    struct mu_to_ell4
      {
        
        double operator()(mu_power n)
          {
            switch(n)
              {
                case mu_power::mu0:
                  return 0.0;
                case mu_power::mu2:
                  return 0.0;
                case mu_power::mu4:
                  return 16.0 / 315.0;
                case mu_power::mu6:
                  return 16.0 / 231.0;
                case mu_power::mu8:
                  return 32.0 / 429.0;
              }
          }
        
      };
    
    struct mu_to_ell4_rs
      {
        
        mu_to_ell4_rs(double _A, double _B)
          : A(_A),
            B(_B)
          {
          }
        
        double operator()(mu_power n)
          {
            switch(n)
              {
                case mu_power::mu0:
                  return (std::exp(-A - B) * (-10. * std::sqrt(A) * (21. + 2. * A) +
                                              3. * (35. + 4. * (-5. + A) * A) * std::exp(A) * std::sqrt(M_PI) *
                                              std::erf(std::sqrt(A)))) / (32. * std::pow(A, 2.5));
                case mu_power::mu2:
                  return (std::exp(-A - B) * (-2. * std::sqrt(A) * (525. + 2. * A * (85. + 16. * A)) +
                                              3. * (175. + 4. * (-15. + A) * A) * std::exp(A) * std::sqrt(M_PI) *
                                              std::erf(std::sqrt(A)))) / (64. * std::pow(A, 3.5));
                case mu_power::mu4:
                  return (std::exp(-A - B) * (-2. * std::sqrt(A) * (3675. + 2. * A * (775. + 16. * A * (13. + 2. * A))) +
                                              3. * (1225. + 12. * (-25. + A) * A) * std::exp(A) * std::sqrt(M_PI) *
                                              std::erf(std::sqrt(A)))) / (128. * std::pow(A, 4.5));
                case mu_power::mu6:
                  return (std::exp(-A - B) *
                          (-2. * std::sqrt(A) * (33075. + 2. * A * (7875. + 32. * A * (75. + A * (15. + 2. * A)))) +
                           45. * (735. + 4. * (-35. + A) * A) * std::exp(A) * std::sqrt(M_PI) * std::erf(std::sqrt(A)))) /
                         (256. * std::pow(A, 5.5));
                case mu_power::mu8:
                  return (std::exp(-A - B) * (-2. * std::sqrt(A) * (363825. + 2. * A * (92925. + 32. * A * (945. + 2. * A * (105. + A * (17. + 2. * A))))) +
                                              315. * (1155. + 4. * (-45. + A) * A) * std::exp(A) * std::sqrt(M_PI) *
                                              std::erf(std::sqrt(A)))) / (512. * std::pow(A, 6.5));
              }
          }
        
        double A;
        double B;
        
      };
    
    
    static int matsubara_A_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        multipole_Pk_calculator_impl::integrand_data* data = static_cast<multipole_Pk_calculator_impl::integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        
        f[0] = data->jacobian * data->k*data->k * data->Pk(q);
        
        return(0);  // return value irrelevant unless = -999, which means stop integration
      }
    
  }

multipole_Pk multipole_Pk_calculator::calculate(const Mpc_units::energy& k, const Mpc_units::energy& IR_resum,
                                                const IR_resum_token& IR_resum_tok, const oneloop_Pk& data,
                                                const oneloop_growth_record& gf_data, const tree_power_spectrum& Ptree)
  {
    // construct lambdas to access components of an RSD P(k) record

    auto access_tree   = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_tree().value; };
    auto access_13     = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_13().value; };
    auto access_22     = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_22().value; };
    auto access_SPT    = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_1loop_SPT().value; };
    auto access_Z2d    = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_delta().value; };
    auto access_Z0v    = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_Z0_v().value; };
    auto access_Z2v    = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_v().value; };
    auto access_Z0vd   = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_Z0_vdelta().value; };
    auto access_Z2vd   = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_vdelta().value; };
    auto access_Z2vv   = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_vv().value; };
    auto access_Z2vvd  = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_vvdelta().value; };
    auto access_Z2vvv  = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_vvv().value; };
    
    // get Matsubara A-coeff (remember we compute it multiplied by k^2 to get something dimensionless)
    double Matsubara_A = this->compute_Matsubara_A(k, IR_resum, Ptree);
    
    double A_coeff = gf_data.f*(gf_data.f+2.0) * Matsubara_A;
    double B_coeff = Matsubara_A;
    
    // compute un-resummed multipole power spectra
    
    Mpc_units::inverse_energy3 P0_tree = this->decompose(access_tree, data, multipole_Pk_calculator_impl::mu_to_ell0());
    Mpc_units::inverse_energy3 P2_tree = this->decompose(access_tree, data, multipole_Pk_calculator_impl::mu_to_ell2());
    Mpc_units::inverse_energy3 P4_tree = this->decompose(access_tree, data, multipole_Pk_calculator_impl::mu_to_ell4());
    
    Mpc_units::inverse_energy3 P0_13 = this->decompose(access_13, data, multipole_Pk_calculator_impl::mu_to_ell0());
    Mpc_units::inverse_energy3 P2_13 = this->decompose(access_13, data, multipole_Pk_calculator_impl::mu_to_ell2());
    Mpc_units::inverse_energy3 P4_13 = this->decompose(access_13, data, multipole_Pk_calculator_impl::mu_to_ell4());
    
    Mpc_units::inverse_energy3 P0_22 = this->decompose(access_22, data, multipole_Pk_calculator_impl::mu_to_ell0());
    Mpc_units::inverse_energy3 P2_22 = this->decompose(access_22, data, multipole_Pk_calculator_impl::mu_to_ell2());
    Mpc_units::inverse_energy3 P4_22 = this->decompose(access_22, data, multipole_Pk_calculator_impl::mu_to_ell4());
    
    Mpc_units::inverse_energy3 P0_SPT = this->decompose(access_SPT, data, multipole_Pk_calculator_impl::mu_to_ell0());
    Mpc_units::inverse_energy3 P2_SPT = this->decompose(access_SPT, data, multipole_Pk_calculator_impl::mu_to_ell2());
    Mpc_units::inverse_energy3 P4_SPT = this->decompose(access_SPT, data, multipole_Pk_calculator_impl::mu_to_ell4());

    // compute resummed multipole power spectra
    
    Mpc_units::inverse_energy3 P0_tree_rs = this->decompose(access_tree, data, multipole_Pk_calculator_impl::mu_to_ell0_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P2_tree_rs = this->decompose(access_tree, data, multipole_Pk_calculator_impl::mu_to_ell2_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P4_tree_rs = this->decompose(access_tree, data, multipole_Pk_calculator_impl::mu_to_ell4_rs(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy3 P0_13_rs = this->decompose(access_13, data, multipole_Pk_calculator_impl::mu_to_ell0_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P2_13_rs = this->decompose(access_13, data, multipole_Pk_calculator_impl::mu_to_ell2_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P4_13_rs = this->decompose(access_13, data, multipole_Pk_calculator_impl::mu_to_ell4_rs(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy3 P0_22_rs = this->decompose(access_22, data, multipole_Pk_calculator_impl::mu_to_ell0_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P2_22_rs = this->decompose(access_22, data, multipole_Pk_calculator_impl::mu_to_ell2_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P4_22_rs = this->decompose(access_22, data, multipole_Pk_calculator_impl::mu_to_ell4_rs(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy3 P0_SPT_rs = this->decompose_1loop_resummed(data, Matsubara_A, gf_data.f, multipole_Pk_calculator_impl::mu_to_ell0_rs(A_coeff, B_coeff), k, Ptree);
    Mpc_units::inverse_energy3 P2_SPT_rs = this->decompose_1loop_resummed(data, Matsubara_A, gf_data.f, multipole_Pk_calculator_impl::mu_to_ell2_rs(A_coeff, B_coeff), k, Ptree);
    Mpc_units::inverse_energy3 P4_SPT_rs = this->decompose_1loop_resummed(data, Matsubara_A, gf_data.f, multipole_Pk_calculator_impl::mu_to_ell4_rs(A_coeff, B_coeff), k, Ptree);
    
    Mpc_units::inverse_energy P0_Z2d_rs = this->decompose(access_Z2d, data, multipole_Pk_calculator_impl::mu_to_ell0_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy P2_Z2d_rs = this->decompose(access_Z2d, data, multipole_Pk_calculator_impl::mu_to_ell2_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy P4_Z2d_rs = this->decompose(access_Z2d, data, multipole_Pk_calculator_impl::mu_to_ell4_rs(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy3 P0_Z0v_rs = this->decompose(access_Z0v, data, multipole_Pk_calculator_impl::mu_to_ell0_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P2_Z0v_rs = this->decompose(access_Z0v, data, multipole_Pk_calculator_impl::mu_to_ell2_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P4_Z0v_rs = this->decompose(access_Z0v, data, multipole_Pk_calculator_impl::mu_to_ell4_rs(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy P0_Z2v_rs = this->decompose(access_Z2v, data, multipole_Pk_calculator_impl::mu_to_ell0_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy P2_Z2v_rs = this->decompose(access_Z2v, data, multipole_Pk_calculator_impl::mu_to_ell2_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy P4_Z2v_rs = this->decompose(access_Z2v, data, multipole_Pk_calculator_impl::mu_to_ell4_rs(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy3 P0_Z0vd_rs = this->decompose(access_Z0vd, data, multipole_Pk_calculator_impl::mu_to_ell0_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P2_Z0vd_rs = this->decompose(access_Z0vd, data, multipole_Pk_calculator_impl::mu_to_ell2_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P4_Z0vd_rs = this->decompose(access_Z0vd, data, multipole_Pk_calculator_impl::mu_to_ell4_rs(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy P0_Z2vd_rs = this->decompose(access_Z2vd, data, multipole_Pk_calculator_impl::mu_to_ell0_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy P2_Z2vd_rs = this->decompose(access_Z2vd, data, multipole_Pk_calculator_impl::mu_to_ell2_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy P4_Z2vd_rs = this->decompose(access_Z2vd, data, multipole_Pk_calculator_impl::mu_to_ell4_rs(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy P0_Z2vv_rs = this->decompose(access_Z2vv, data, multipole_Pk_calculator_impl::mu_to_ell0_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy P2_Z2vv_rs = this->decompose(access_Z2vv, data, multipole_Pk_calculator_impl::mu_to_ell2_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy P4_Z2vv_rs = this->decompose(access_Z2vv, data, multipole_Pk_calculator_impl::mu_to_ell4_rs(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy P0_Z2vvd_rs = this->decompose(access_Z2vvd, data, multipole_Pk_calculator_impl::mu_to_ell0_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy P2_Z2vvd_rs = this->decompose(access_Z2vvd, data, multipole_Pk_calculator_impl::mu_to_ell2_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy P4_Z2vvd_rs = this->decompose(access_Z2vvd, data, multipole_Pk_calculator_impl::mu_to_ell4_rs(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy P0_Z2vvv_rs = this->decompose(access_Z2vvv, data, multipole_Pk_calculator_impl::mu_to_ell0_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy P2_Z2vvv_rs = this->decompose(access_Z2vvv, data, multipole_Pk_calculator_impl::mu_to_ell2_rs(A_coeff, B_coeff));
    Mpc_units::inverse_energy P4_Z2vvv_rs = this->decompose(access_Z2vvv, data, multipole_Pk_calculator_impl::mu_to_ell4_rs(A_coeff, B_coeff));
    
    Pk_ell P0(P0_tree, P0_tree_rs, P0_13, P0_13_rs, P0_22, P0_22_rs, P0_SPT, P0_SPT_rs,
              P0_Z2d_rs, P0_Z0v_rs, P0_Z2v_rs, P0_Z0vd_rs, P0_Z2vd_rs, P0_Z2vv_rs, P0_Z2vvd_rs, P0_Z2vvv_rs);

    Pk_ell P2(P2_tree, P2_tree_rs, P2_13, P2_13_rs, P2_22, P2_22_rs, P2_SPT, P2_SPT_rs,
              P2_Z2d_rs, P2_Z0v_rs, P2_Z2v_rs, P2_Z0vd_rs, P2_Z2vd_rs, P2_Z2vv_rs, P2_Z2vvd_rs, P2_Z2vvv_rs);

    Pk_ell P4(P4_tree, P4_tree_rs, P4_13, P4_13_rs, P4_22, P4_22_rs, P4_SPT, P4_SPT_rs,
              P4_Z2d_rs, P4_Z0v_rs, P4_Z2v_rs, P4_Z0vd_rs, P4_Z2vd_rs, P4_Z2vv_rs, P4_Z2vvd_rs, P4_Z2vvv_rs);
    
    return multipole_Pk(data.get_k_token(), data.get_IR_token(),
                        data.get_UV_token(), data.get_z_token(),
                        IR_resum_tok, P0, P2, P4);
  }


template <typename Accessor, typename Decomposer>
decltype(std::declval<Accessor>()(std::declval<const rsd_dd_Pk&>()))
multipole_Pk_calculator::decompose(Accessor extract, const oneloop_Pk& data, Decomposer decomp)
  {
    auto mu0 = extract(data.get_dd_rsd_mu0());
    auto mu2 = extract(data.get_dd_rsd_mu2());
    auto mu4 = extract(data.get_dd_rsd_mu4());
    auto mu6 = extract(data.get_dd_rsd_mu6());
    auto mu8 = extract(data.get_dd_rsd_mu8());
    
    return mu0 * decomp(mu_power::mu0)
           + mu2 * decomp(mu_power::mu2)
           + mu4 * decomp(mu_power::mu4)
           + mu6 * decomp(mu_power::mu6)
           + mu8 * decomp(mu_power::mu8);
  }


template <typename Decomposer>
Mpc_units::inverse_energy3
multipole_Pk_calculator::decompose_1loop_resummed(const oneloop_Pk& data, double Matsubara_A, double f,
                                                  Decomposer decomp, const Mpc_units::energy& k,
                                                  const tree_power_spectrum& Ptree)
  {
    auto raw_mu0 = data.get_dd_rsd_mu0().get_1loop_SPT().value;
    auto raw_mu2 = data.get_dd_rsd_mu2().get_1loop_SPT().value;
    auto raw_mu4 = data.get_dd_rsd_mu4().get_1loop_SPT().value;
    auto raw_mu6 = data.get_dd_rsd_mu6().get_1loop_SPT().value;
    auto raw_mu8 = data.get_dd_rsd_mu8().get_1loop_SPT().value;
    
    auto mu0 = raw_mu0 + Matsubara_A * Ptree(k);
    auto mu2 = raw_mu2 + (4.0*f + f*f) * Matsubara_A * Ptree(k);
    auto mu4 = raw_mu4 + (5.0*f*f + 2.0*f*f*f) * Matsubara_A * Ptree(k);
    auto mu6 = raw_mu6 + (2.0*f*f*f + f*f*f*f) * Matsubara_A * Ptree(k);
    auto mu8 = raw_mu8;
    
    return mu0 * decomp(mu_power::mu0)
           + mu2 * decomp(mu_power::mu2)
           + mu4 * decomp(mu_power::mu4)
           + mu6 * decomp(mu_power::mu6)
           + mu8 * decomp(mu_power::mu8);
  }


double multipole_Pk_calculator::compute_Matsubara_A(const Mpc_units::energy& k, const Mpc_units::energy& IR_resum,
                                                    const tree_power_spectrum& Ptree)
  {
    cubareal integral[multipole_Pk_calculator_impl::dimensions];
    cubareal error[multipole_Pk_calculator_impl::dimensions];
    cubareal prob[multipole_Pk_calculator_impl::dimensions];
    
    int regions;
    int evaluations;
    int fail;
    
    const powerspectrum_database& db = Ptree.get_db();
    
    std::unique_ptr<multipole_Pk_calculator_impl::integrand_data> data =
      std::make_unique<multipole_Pk_calculator_impl::integrand_data>(k, IR_resum, db.get_k_min(), Ptree);
    
    cubacores(0, multipole_Pk_calculator_impl::pcores);
    
    Cuhre(multipole_Pk_calculator_impl::dimensions,
          multipole_Pk_calculator_impl::components,
          multipole_Pk_calculator_impl::matsubara_A_integrand, data.get(),
          multipole_Pk_calculator_impl::points_per_invocation,
          this->rel_err, this->abs_err,
          multipole_Pk_calculator_impl::verbosity_none | multipole_Pk_calculator_impl::samples_last,
          multipole_Pk_calculator_impl::min_eval, multipole_Pk_calculator_impl::max_eval,
          multipole_Pk_calculator_impl::cuhre_key,
          nullptr, nullptr,
          &regions, &evaluations, &fail,
          integral, error, prob);
    
    return integral[0] / (6.0*M_PI*M_PI);
  }
