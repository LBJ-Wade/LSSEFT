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
        
        integrand_data(const Mpc_units::energy& UV, const Mpc_units::energy& IR, const tree_power_spectrum& _Pk)
          : UV_cutoff(UV),
            IR_cutoff(IR),
            Pk(_Pk),
            jacobian(UV_cutoff - IR_cutoff),
            q_range(UV_cutoff - IR_cutoff)
          {
          }
        
        const Mpc_units::energy& UV_cutoff;
        const Mpc_units::energy& IR_cutoff;
        const tree_power_spectrum& Pk;
        
        Mpc_units::energy  jacobian;
        Mpc_units::energy  q_range;
      };
    
    
    struct mu_to_ell0
      {
        
        double operator()(mu_power n)
          {
            switch(n)
              {
                case mu_power::mu0: return 2.0;
                case mu_power::mu2: return 2.0 / 3.0;
                case mu_power::mu4: return 2.0 / 5.0;
                case mu_power::mu6: return 2.0 / 7.0;
                case mu_power::mu8: return 2.0 / 9.0;
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
        
        double mu0()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (2 - (2*A)/3. + std::pow(A,2)/5. - std::pow(A,3)/21. + std::pow(A,4)/108. - std::pow(A,5)/660. + std::pow(A,6)/4680. - std::pow(A,7)/37800. + std::pow(A,8)/342720.);
            
            return (std::sqrt(M_PI)*std::erf(std::sqrt(A)))/(std::sqrt(A)*std::exp(B));
          }
        
        double mu2()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.6666666666666666 - (2*A)/5. + std::pow(A,2)/7. - std::pow(A,3)/27. + std::pow(A,4)/132. - std::pow(A,5)/780. + std::pow(A,6)/5400. - std::pow(A,7)/42840. +
                                                          std::pow(A,8)/383040.);
            
            return (std::exp(-A - B)*(-2*std::sqrt(A) + std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(2.*std::pow(A,1.5));
          }
        
        double mu4()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.4 - (2*A)/7. + std::pow(A,2)/9. - std::pow(A,3)/33. + std::pow(A,4)/156. - std::pow(A,5)/900. + std::pow(A,6)/6120. - std::pow(A,7)/47880. +
                                                          std::pow(A,8)/423360.);
            
            return (std::exp(-A - B)*(-2*std::sqrt(A)*(3 + 2*A) + 3*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(4.*std::pow(A,2.5));
          }
        
        double mu6()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.2857142857142857 - (2*A)/9. + std::pow(A,2)/11. - std::pow(A,3)/39. + std::pow(A,4)/180. - std::pow(A,5)/1020. + std::pow(A,6)/6840. - std::pow(A,7)/52920. +
                                                          std::pow(A,8)/463680.);
            
            return (std::exp(-A - B)*(-2*std::sqrt(A)*(15 + 2*A*(5 + 2*A)) + 15*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(8.*std::pow(A,3.5));
          }
        
        double mu8()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.2222222222222222 - (2*A)/11. + std::pow(A,2)/13. - std::pow(A,3)/45. + std::pow(A,4)/204. - std::pow(A,5)/1140. + std::pow(A,6)/7560. -
                                                          std::pow(A,7)/57960. + std::pow(A,8)/504000.);
            
            return (std::exp(-A - B)*(-2*std::sqrt(A)*(105 + 2*A*(35 + 2*A*(7 + 2*A))) + 105*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(16.*std::pow(A,4.5));
          }
        
        double operator()(mu_power n)
          {
            switch(n)
              {
                case mu_power::mu0: return this->mu0();
                case mu_power::mu2: return this->mu2();
                case mu_power::mu4: return this->mu4();
                case mu_power::mu6: return this->mu6();
                case mu_power::mu8: return this->mu8();
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
                case mu_power::mu0: return 0.0;
                case mu_power::mu2: return 4.0 / 15.0;
                case mu_power::mu4: return 8.0 / 35.0;
                case mu_power::mu6: return 4.0 / 21.0;
                case mu_power::mu8: return 16.0 / 99.0;
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
    
        double mu0()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * ((-4*A)/15. + (4*std::pow(A,2))/35. - (2*std::pow(A,3))/63. + (2*std::pow(A,4))/297. - std::pow(A,5)/858. + std::pow(A,6)/5850. - std::pow(A,7)/45900. +
                                                          std::pow(A,8)/406980.);
        
            return (std::exp(-A - B)*(-6*std::sqrt(A) - (-3 + 2*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(4.*std::pow(A,1.5));
          }
    
        double mu2()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.26666666666666666 - (8*A)/35. + (2*std::pow(A,2))/21. - (8*std::pow(A,3))/297. + (5*std::pow(A,4))/858. - std::pow(A,5)/975. + (7*std::pow(A,6))/45900. -
                                                          (2*std::pow(A,7))/101745. + std::pow(A,8)/446880.);
        
            return (std::exp(-A - B)*(-2*std::sqrt(A)*(9 + 4*A) - (-9 + 2*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(8.*std::pow(A,2.5));
          }
    
        double mu4()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.22857142857142856 - (4*A)/21. + (8*std::pow(A,2))/99. - (10*std::pow(A,3))/429. + std::pow(A,4)/195. - (7*std::pow(A,5))/7650. +
                                                          (2*std::pow(A,6))/14535. - std::pow(A,7)/55860. + std::pow(A,8)/486864.);
        
            return (std::exp(-A - B)*(-2*std::sqrt(A)*(45 + 8*A*(3 + A)) - 3*(-15 + 2*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(16.*std::pow(A,3.5));
          }
    
        double mu6()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.19047619047619047 - (16*A)/99. + (10*std::pow(A,2))/143. - (4*std::pow(A,3))/195. + (7*std::pow(A,4))/1530. - (4*std::pow(A,5))/4845. +
                                                          std::pow(A,6)/7980. - std::pow(A,7)/60858. + (11*std::pow(A,8))/5.796e6);
        
            return (std::exp(-A - B)*(-2*std::sqrt(A)*(315 + 4*A*(45 + 4*A*(4 + A))) - 15*(-21 + 2*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(32.*std::pow(A,4.5));
          }
    
        double mu8()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.16161616161616163 - (20*A)/143. + (4*std::pow(A,2))/65. - (14*std::pow(A,3))/765. + (4*std::pow(A,4))/969. - std::pow(A,5)/1330. + std::pow(A,6)/8694. -
                                                          (11*std::pow(A,7))/724500. + std::pow(A,8)/567000.);
        
            return (std::exp(-A - B)*(-2*std::sqrt(A)*(2835 + 8*A*(210 + A*(77 + 4*A*(5 + A)))) - 105*(-27 + 2*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(64.*std::pow(A,5.5));
          }
    
        double operator()(mu_power n)
          {
            switch(n)
              {
                case mu_power::mu0: return this->mu0();
                case mu_power::mu2: return this->mu2();
                case mu_power::mu4: return this->mu4();
                case mu_power::mu6: return this->mu6();
                case mu_power::mu8: return this->mu8();
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
                case mu_power::mu0: return 0.0;
                case mu_power::mu2: return 0.0;
                case mu_power::mu4: return 16.0 / 315.0;
                case mu_power::mu6: return 16.0 / 231.0;
                case mu_power::mu8: return 32.0 / 429.0;
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
    
        double mu0()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * ((8*std::pow(A,2))/315. - (8*std::pow(A,3))/693. + (4*std::pow(A,4))/1287. - (4*std::pow(A,5))/6435. + std::pow(A,6)/9945. - std::pow(A,7)/72675. +
                                                          std::pow(A,8)/610470.);
        
            return (std::exp(-A - B)*(-10*std::sqrt(A)*(21 + 2*A) + 3*(35 + 4*(-5 + A)*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(32.*std::pow(A,2.5));
          }
    
        double mu2()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * ((-16*A)/315. + (8*std::pow(A,2))/231. - (16*std::pow(A,3))/1287. + (4*std::pow(A,4))/1287. - (2*std::pow(A,5))/3315. + (7*std::pow(A,6))/72675. -
                                                          (4*std::pow(A,7))/305235. + std::pow(A,8)/642390.);
        
            return (std::exp(-A - B)*(-2*std::sqrt(A)*(525 + 2*A*(85 + 16*A)) + 3*(175 + 4*(-15 + A)*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(64.*std::pow(A,3.5));
          }
    
        double mu4()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.050793650793650794 - (16*A)/231. + (16*std::pow(A,2))/429. - (16*std::pow(A,3))/1287. + (2*std::pow(A,4))/663. - (14*std::pow(A,5))/24225. +
                                                          (4*std::pow(A,6))/43605. - (4*std::pow(A,7))/321195. + std::pow(A,8)/676200.);
        
            return (std::exp(-A - B)*(-2*std::sqrt(A)*(3675 + 2*A*(775 + 16*A*(13 + 2*A))) + 3*(1225 + 12*(-25 + A)*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/
                   (128.*std::pow(A,4.5));
          }
    
        double mu6()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.06926406926406926 - (32*A)/429. + (16*std::pow(A,2))/429. - (8*std::pow(A,3))/663. + (14*std::pow(A,4))/4845. - (8*std::pow(A,5))/14535. +
                                                          (4*std::pow(A,6))/45885. - std::pow(A,7)/84525. + (11*std::pow(A,8))/7.8246e6);
        
            return (std::exp(-A - B)*(-2*std::sqrt(A)*(33075 + 2*A*(7875 + 32*A*(75 + A*(15 + 2*A)))) + 45*(735 + 4*(-35 + A)*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/
                   (256.*std::pow(A,5.5));
          }
    
        double mu8()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.07459207459207459 - (32*A)/429. + (8*std::pow(A,2))/221. - (56*std::pow(A,3))/4845. + (8*std::pow(A,4))/2907. - (8*std::pow(A,5))/15295. +
                                                          std::pow(A,6)/12075. - (11*std::pow(A,7))/978075. + (11*std::pow(A,8))/8.2215e6);
        
            return (std::exp(-A - B)*(-2*std::sqrt(A)*(363825 + 2*A*(92925 + 32*A*(945 + 2*A*(105 + A*(17 + 2*A))))) +
                                     315*(1155 + 4*(-45 + A)*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(512.*std::pow(A,6.5));
          }
    
        double operator()(mu_power n)
          {
            switch(n)
              {
                case mu_power::mu0: return this->mu0();
                case mu_power::mu2: return this->mu2();
                case mu_power::mu4: return this->mu4();
                case mu_power::mu6: return this->mu6();
                case mu_power::mu8: return this->mu8();
              }
          }
        
        double A;
        double B;
        
      };
    
    
    static int matsubara_A_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        multipole_Pk_calculator_impl::integrand_data* data = static_cast<multipole_Pk_calculator_impl::integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        
        f[0] = data->jacobian * data->Pk(q) / Mpc_units::Mpc2;
        
        return(0);  // return value irrelevant unless = -999, which means stop integration
      }
    
  }

multipole_Pk multipole_Pk_calculator::calculate_Legendre(const Mpc_units::energy& k, const Matsubara_A& A,
                                                         const oneloop_Pk& data, const oneloop_growth_record& gf_data,
                                                         const tree_power_spectrum& Ptree)
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
    
    // get Matsubara A-coeff (remember we have to scale up by the linear growth factor, since we store just the
    // raw integral over the early-time tree-level power spectrum)
    double Matsubara_A = gf_data.g*gf_data.g * k*k * A;
    
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
    
    Mpc_units::inverse_energy3 P0_SPT_rs = this->decompose_1loop_resummed(data, Matsubara_A, gf_data, multipole_Pk_calculator_impl::mu_to_ell0_rs(A_coeff, B_coeff), k, Ptree);
    Mpc_units::inverse_energy3 P2_SPT_rs = this->decompose_1loop_resummed(data, Matsubara_A, gf_data, multipole_Pk_calculator_impl::mu_to_ell2_rs(A_coeff, B_coeff), k, Ptree);
    Mpc_units::inverse_energy3 P4_SPT_rs = this->decompose_1loop_resummed(data, Matsubara_A, gf_data, multipole_Pk_calculator_impl::mu_to_ell4_rs(A_coeff, B_coeff), k, Ptree);
    
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
                        A.get_token(), P0, P2, P4);
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
multipole_Pk_calculator::decompose_1loop_resummed(const oneloop_Pk& data, double Matsubara_A,
                                                  const oneloop_growth_record& gf_data,
                                                  Decomposer decomp, const Mpc_units::energy& k,
                                                  const tree_power_spectrum& Ptree)
  {
    auto raw_mu0 = data.get_dd_rsd_mu0().get_1loop_SPT().value;
    auto raw_mu2 = data.get_dd_rsd_mu2().get_1loop_SPT().value;
    auto raw_mu4 = data.get_dd_rsd_mu4().get_1loop_SPT().value;
    auto raw_mu6 = data.get_dd_rsd_mu6().get_1loop_SPT().value;
    auto raw_mu8 = data.get_dd_rsd_mu8().get_1loop_SPT().value;
    
    double g = gf_data.g;
    double f = gf_data.f;
    
    auto Ptr = g*g * Ptree(k);
    
    auto mu0 = raw_mu0 + Matsubara_A * Ptr;
    auto mu2 = raw_mu2 + (4.0*f + f*f) * Matsubara_A * Ptr;
    auto mu4 = raw_mu4 + (5.0*f*f + 2.0*f*f*f) * Matsubara_A * Ptr;
    auto mu6 = raw_mu6 + (2.0*f*f*f + f*f*f*f) * Matsubara_A * Ptr;
    auto mu8 = raw_mu8;
    
    return mu0 * decomp(mu_power::mu0)
           + mu2 * decomp(mu_power::mu2)
           + mu4 * decomp(mu_power::mu4)
           + mu6 * decomp(mu_power::mu6)
           + mu8 * decomp(mu_power::mu8);
  }


Matsubara_A
multipole_Pk_calculator::calculate_Matsubara_A(const Mpc_units::energy& IR_resum, const IR_resum_token& IR_resum_tok,
                                               const tree_power_spectrum& Ptree)
  {
    cubareal integral[multipole_Pk_calculator_impl::dimensions];
    cubareal error[multipole_Pk_calculator_impl::dimensions];
    cubareal prob[multipole_Pk_calculator_impl::dimensions];
    
    int regions;
    int evaluations;
    int fail;
    
    const powerspectrum_database& db = Ptree.get_db();
    
    // use 10% clearance above lower limit of spline to avoid unwanted effects associated
    // with inaccuracies in the fit there
    constexpr double TEN_PERCENT_CLEARANCE = 1.1;
    
    std::unique_ptr<multipole_Pk_calculator_impl::integrand_data> data =
      std::make_unique<multipole_Pk_calculator_impl::integrand_data>(IR_resum, TEN_PERCENT_CLEARANCE*db.get_k_min(), Ptree);
    
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
    
    return Matsubara_A(IR_resum_tok, integral[0] * Mpc_units::Mpc2 / (6.0*M_PI*M_PI));
  }
