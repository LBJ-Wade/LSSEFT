//
// Created by David Seery on 18/11/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#include <cmath>

#include "cosmology/concepts/power_spectrum.h"

#include "multipole_Pk_calculator.h"

#include "cuba.h"
#include "boost/math/special_functions/bessel.hpp"

namespace multipole_Pk_calculator_impl
  {

    constexpr unsigned int dimensions            = 2;       // k and q integrals
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
        
        integrand_data(const Mpc_units::energy& UV, const Mpc_units::energy& IR,
                       const Mpc_units::inverse_energy& _qmin, const Mpc_units::inverse_energy& _qmax,
                       const spline_Pk& _Pk)
          : UV_cutoff(UV),
            IR_cutoff(IR),
            qmin(_qmin),
            qmax(_qmax),
            Pk(_Pk),
            jacobian((UV_cutoff - IR_cutoff) * (qmax - qmin)),
            s_range(UV_cutoff - IR_cutoff),
            q_range(qmax - qmin),
            q_volume(qmax*qmax*qmax/3.0 - qmin*qmin*qmin/3.0)
          {
          }
        
        const Mpc_units::energy& UV_cutoff;
        const Mpc_units::energy& IR_cutoff;
        
        const Mpc_units::inverse_energy& qmin;
        const Mpc_units::inverse_energy& qmax;
        
        const spline_Pk& Pk;
        
        double jacobian;

        Mpc_units::energy s_range;
        Mpc_units::inverse_energy q_range;
        Mpc_units::inverse_energy3 q_volume;
      };
    
    
    struct mu_to_ell0
      {
        
        double operator()(mu_power n)
          {
            switch(n)
              {
                case mu_power::mu0: return 1.0;
                case mu_power::mu2: return 1.0 / 3.0;
                case mu_power::mu4: return 1.0 / 5.0;
                case mu_power::mu6: return 1.0 / 7.0;
                case mu_power::mu8: return 1.0 / 9.0;
              }
          }
        
      };
    
    struct mu_to_ell0_expXY
      {
        
        mu_to_ell0_expXY(double _A, double _B)
          : A(_A),
            B(_B)
          {
          }
        
        double mu0()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (1 - A/3. + std::pow(A,2)/10. - std::pow(A,3)/42. + std::pow(A,4)/216. - std::pow(A,5)/1320. + std::pow(A,6)/9360. - std::pow(A,7)/75600. + std::pow(A,8)/685440.);
            
            return (std::sqrt(M_PI)*std::erf(std::sqrt(A)))/(2.*std::sqrt(A)*std::exp(B));
          }
        
        double mu2()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.3333333333333333 - A/5. + std::pow(A,2)/14. - std::pow(A,3)/54. + std::pow(A,4)/264. - std::pow(A,5)/1560. + std::pow(A,6)/10800. - std::pow(A,7)/85680. +
                                                          std::pow(A,8)/766080.);
            
            return (std::exp(-A - B)*(-2*std::sqrt(A) + std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(4.*std::pow(A,1.5));
          }
        
        double mu4()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.2 - A/7. + std::pow(A,2)/18. - std::pow(A,3)/66. + std::pow(A,4)/312. - std::pow(A,5)/1800. + std::pow(A,6)/12240. - std::pow(A,7)/95760. + std::pow(A,8)/846720.);
            
            return (std::exp(-A - B)*(-2*std::sqrt(A)*(3 + 2*A) + 3*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(8.*std::pow(A,2.5));
          }
        
        double mu6()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.14285714285714285 - A/9. + std::pow(A,2)/22. - std::pow(A,3)/78. + std::pow(A,4)/360. - std::pow(A,5)/2040. + std::pow(A,6)/13680. - std::pow(A,7)/105840. +
                                                          std::pow(A,8)/927360.);
            
            return (std::exp(-A - B)*(-2*std::sqrt(A)*(15 + 2*A*(5 + 2*A)) + 15*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(16.*std::pow(A,3.5));
          }
        
        double mu8()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.1111111111111111 - A/11. + std::pow(A,2)/26. - std::pow(A,3)/90. + std::pow(A,4)/408. - std::pow(A,5)/2280. + std::pow(A,6)/15120. - std::pow(A,7)/115920. +
                                                          std::pow(A,8)/1.008e6);
            
            return (std::exp(-A - B)*(-2*std::sqrt(A)*(105 + 2*A*(35 + 2*A*(7 + 2*A))) + 105*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(32.*std::pow(A,4.5));
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
                case mu_power::mu2: return 2.0 / 3.0;
                case mu_power::mu4: return 4.0 / 7.0;
                case mu_power::mu6: return 10.0 / 21.0;
                case mu_power::mu8: return 40.0 / 99.0;
              }
          }
        
      };
    
    struct mu_to_ell2_expXY
      {
        
        mu_to_ell2_expXY(double _A, double _B)
          : A(_A),
            B(_B)
          {
          }
    
        double mu0()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * ((-2*A)/3. + (2*std::pow(A,2))/7. - (5*std::pow(A,3))/63. + (5*std::pow(A,4))/297. - (5*std::pow(A,5))/1716. + std::pow(A,6)/2340. - std::pow(A,7)/18360. +
                                                          std::pow(A,8)/162792.);
        
            return (-5*std::exp(-A - B)*(6*std::sqrt(A) + (-3 + 2*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(8.*std::pow(A,1.5));
          }
    
        double mu2()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.6666666666666666 - (4*A)/7. + (5*std::pow(A,2))/21. - (20*std::pow(A,3))/297. + (25*std::pow(A,4))/1716. - std::pow(A,5)/390. + (7*std::pow(A,6))/18360. -
                                                          std::pow(A,7)/20349. + std::pow(A,8)/178752.);
        
            return (-5*std::exp(-A - B)*(2*std::sqrt(A)*(9 + 4*A) + (-9 + 2*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(16.*std::pow(A,2.5));
          }
    
        double mu4()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.5714285714285714 - (10*A)/21. + (20*std::pow(A,2))/99. - (25*std::pow(A,3))/429. + std::pow(A,4)/78. - (7*std::pow(A,5))/3060. + std::pow(A,6)/2907. -
                                                          std::pow(A,7)/22344. + (5*std::pow(A,8))/973728.);
        
            return (-5*std::exp(-A - B)*(2*std::sqrt(A)*(45 + 8*A*(3 + A)) + 3*(-15 + 2*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(32.*std::pow(A,3.5));
          }
    
        double mu6()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.47619047619047616 - (40*A)/99. + (25*std::pow(A,2))/143. - (2*std::pow(A,3))/39. + (7*std::pow(A,4))/612. - (2*std::pow(A,5))/969. + std::pow(A,6)/3192. -
                                                          (5*std::pow(A,7))/121716. + (11*std::pow(A,8))/2.3184e6);
        
            return (std::exp(-A - B)*(-10*std::sqrt(A)*(315 + 4*A*(45 + 4*A*(4 + A))) - 75*(-21 + 2*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(64.*std::pow(A,4.5));
          }
    
        double mu8()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.40404040404040403 - (50*A)/143. + (2*std::pow(A,2))/13. - (7*std::pow(A,3))/153. + (10*std::pow(A,4))/969. - std::pow(A,5)/532. +
                                                          (5*std::pow(A,6))/17388. - (11*std::pow(A,7))/289800. + std::pow(A,8)/226800.);
        
            return (5*std::exp(-A - B)*(-2*std::sqrt(A)*(2835 + 8*A*(210 + A*(77 + 4*A*(5 + A)))) - 105*(-27 + 2*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/
                   (128.*std::pow(A,5.5));
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
                case mu_power::mu4: return 8.0 / 35.0;
                case mu_power::mu6: return 24.0 / 77.0;
                case mu_power::mu8: return 48.0 / 143.0;
              }
          }
        
      };
    
    struct mu_to_ell4_expXY
      {
        
        mu_to_ell4_expXY(double _A, double _B)
          : A(_A),
            B(_B)
          {
          }
    
        double mu0()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * ((4*std::pow(A,2))/35. - (4*std::pow(A,3))/77. + (2*std::pow(A,4))/143. - (2*std::pow(A,5))/715. + std::pow(A,6)/2210. - std::pow(A,7)/16150. +
                                                          std::pow(A,8)/135660.);
        
            return (std::exp(-A - B)*(-90*std::sqrt(A)*(21 + 2*A) + 27*(35 + 4*(-5 + A)*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(64.*std::pow(A,2.5));
          }
    
        double mu2()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * ((-8*A)/35. + (12*std::pow(A,2))/77. - (8*std::pow(A,3))/143. + (2*std::pow(A,4))/143. - (3*std::pow(A,5))/1105. + (7*std::pow(A,6))/16150. -
                                                          (2*std::pow(A,7))/33915. + (3*std::pow(A,8))/428260.);
        
            return (std::exp(-A - B)*(-18*std::sqrt(A)*(525 + 2*A*(85 + 16*A)) + 27*(175 + 4*(-15 + A)*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(128.*std::pow(A,3.5));
          }
    
        double mu4()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.22857142857142856 - (24*A)/77. + (24*std::pow(A,2))/143. - (8*std::pow(A,3))/143. + (3*std::pow(A,4))/221. - (21*std::pow(A,5))/8075. +
                                                          (2*std::pow(A,6))/4845. - (6*std::pow(A,7))/107065. + (3*std::pow(A,8))/450800.);
        
            return (std::exp(-A - B)*(-18*std::sqrt(A)*(3675 + 2*A*(775 + 16*A*(13 + 2*A))) + 27*(1225 + 12*(-25 + A)*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/
                   (256.*std::pow(A,4.5));
          }
    
        double mu6()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.3116883116883117 - (48*A)/143. + (24*std::pow(A,2))/143. - (12*std::pow(A,3))/221. + (21*std::pow(A,4))/1615. - (4*std::pow(A,5))/1615. +
                                                          (6*std::pow(A,6))/15295. - (3*std::pow(A,7))/56350. + (11*std::pow(A,8))/1.7388e6);
        
            return (9*std::exp(-A - B)*(-2*std::sqrt(A)*(33075 + 2*A*(7875 + 32*A*(75 + A*(15 + 2*A)))) + 45*(735 + 4*(-35 + A)*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/
                   (512.*std::pow(A,5.5));
          }
    
        double mu8()
          {
            if(std::abs(A) < 0.01) return std::exp(-B) * (0.3356643356643357 - (48*A)/143. + (36*std::pow(A,2))/221. - (84*std::pow(A,3))/1615. + (4*std::pow(A,4))/323. - (36*std::pow(A,5))/15295. +
                                                          (3*std::pow(A,6))/8050. - (11*std::pow(A,7))/217350. + (11*std::pow(A,8))/1.827e6);
        
            return (9*std::exp(-A - B)*(-2*std::sqrt(A)*(363825 + 2*A*(92925 + 32*A*(945 + 2*A*(105 + A*(17 + 2*A))))) +
                                       315*(1155 + 4*(-45 + A)*A)*std::exp(A)*std::sqrt(M_PI)*std::erf(std::sqrt(A))))/(1024.*std::pow(A,6.5));
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
    
    
    // subtractions for resummation of full 1-loop power spectrum
    struct resum_adjuster
      {

        resum_adjuster(const Mpc_units::energy& _k, double _XY, const oneloop_growth_record& _gf, const wiggle_Pk& _Pk)
          : k(_k),
            XY(_XY),
            gf(_gf),
            Pk(_Pk),
            f(gf.f),
            // the factor appearing in each subtraction is k^2 (X+Y) P_lin,nw
            // here, k^2 (X+Y) is our variable XY
            // P_lin,w is g^2 Pk.Pk_wiggle
            factor(XY * gf.g * gf.g * Pk.Pk_wiggle(k))
          {
          }
    
        Mpc_units::inverse_energy3 operator()(mu_power n)
          {
            switch(n)
              {
                case mu_power::mu0: return factor;
                case mu_power::mu2: return (4.0*f + f*f) * factor;
                case mu_power::mu4: return (5.0*f*f + 2.0*f*f*f) * factor;
                case mu_power::mu6: return (2.0*f*f*f + f*f*f*f) * factor;
                case mu_power::mu8: return 0.0;
              }
          }
        
        const Mpc_units::energy& k;
        double XY;
        const oneloop_growth_record& gf;
        const wiggle_Pk& Pk;
        
        double f;
        Mpc_units::inverse_energy3 factor;
        
      };
    
    
    // null subtractions
    template <typename Dimension>
    struct null_adjuster
      {

        Dimension operator()(mu_power n)
          {
            switch(n)
              {
                case mu_power::mu0: return 0.0;
                case mu_power::mu2: return 0.0;
                case mu_power::mu4: return 0.0;
                case mu_power::mu6: return 0.0;
                case mu_power::mu8: return 0.0;
              }
          }
        
      };
    
    
    static int matsubara_X_integrand(const int* ndim, const cubareal* x, const int* ncomp, cubareal* f, void* userdata)
      {
        multipole_Pk_calculator_impl::integrand_data* data = static_cast<multipole_Pk_calculator_impl::integrand_data*>(userdata);
        
        Mpc_units::energy s = data->IR_cutoff + x[0] * data->s_range;
        Mpc_units::inverse_energy q = data->qmin + x[1] * data->q_range;
        
        Mpc_units::energy qfactor = q*q / data->q_volume;
        
        f[0] = data->jacobian * qfactor * data->Pk(s) * (1.0 - (3.0/(q*s))*boost::math::sph_bessel(1, q*s)) / Mpc_units::Mpc2;
        
        return(0);  // return value irrelevant unless = -999, which means stop integration
      }
    
    
    static int matsubara_Y_integrand(const int* ndim, const cubareal* x, const int* ncomp, cubareal* f, void* userdata)
      {
        multipole_Pk_calculator_impl::integrand_data* data = static_cast<multipole_Pk_calculator_impl::integrand_data*>(userdata);
        
        Mpc_units::energy s = data->IR_cutoff + x[0] * data->s_range;
        Mpc_units::inverse_energy q = data->qmin + x[1] * data->q_range;
    
        Mpc_units::energy qfactor = q*q / data->q_volume;
    
        f[0] = data->jacobian * qfactor * data->Pk(s) * boost::math::sph_bessel(2, q*s) / Mpc_units::Mpc2;
        
        return(0);  // return value irrelevant unless = -999, which means stop integration
      }
    
  }


template <typename Accessor, typename Decomposer>
auto
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


template <typename WiggleAccessor, typename NoWiggleAccessor, typename ResumAdjuster, typename RawDecomposer, typename XYDecomposer>
auto
multipole_Pk_calculator::decompose(WiggleAccessor wiggle, NoWiggleAccessor nowiggle, const oneloop_Pk& data,
                                   ResumAdjuster adjust, RawDecomposer raw_decomp, XYDecomposer XY_decomp)
  {
    // the decomposition consists of two parts
    // the broadband part comes from the no-wiggle power spectrum and has no resummation
    auto nowiggle_mu0 = nowiggle(data.get_dd_rsd_mu0());
    auto nowiggle_mu2 = nowiggle(data.get_dd_rsd_mu2());
    auto nowiggle_mu4 = nowiggle(data.get_dd_rsd_mu4());
    auto nowiggle_mu6 = nowiggle(data.get_dd_rsd_mu6());
    auto nowiggle_mu8 = nowiggle(data.get_dd_rsd_mu8());
    
    auto nowiggle_Pl = nowiggle_mu0 * raw_decomp(mu_power::mu0)
                       + nowiggle_mu2 * raw_decomp(mu_power::mu2)
                       + nowiggle_mu4 * raw_decomp(mu_power::mu4)
                       + nowiggle_mu6 * raw_decomp(mu_power::mu6)
                       + nowiggle_mu8 * raw_decomp(mu_power::mu8);
    
    // meanwhile, the wiggle part is resummed
    // its mu coefficients must be adjusted to account for subtractions associated with the resummation
    auto wiggle_mu0 = wiggle(data.get_dd_rsd_mu0()) + adjust(mu_power::mu0);
    auto wiggle_mu2 = wiggle(data.get_dd_rsd_mu2()) + adjust(mu_power::mu2);
    auto wiggle_mu4 = wiggle(data.get_dd_rsd_mu4()) + adjust(mu_power::mu4);
    auto wiggle_mu6 = wiggle(data.get_dd_rsd_mu6()) + adjust(mu_power::mu6);
    auto wiggle_mu8 = wiggle(data.get_dd_rsd_mu8()) + adjust(mu_power::mu8);
    
    auto wiggle_Pl = wiggle_mu0 * XY_decomp(mu_power::mu0)
                     + wiggle_mu2 * XY_decomp(mu_power::mu2)
                     + wiggle_mu4 * XY_decomp(mu_power::mu4)
                     + wiggle_mu6 * XY_decomp(mu_power::mu6)
                     + wiggle_mu8 * XY_decomp(mu_power::mu8);
    
    return nowiggle_Pl + wiggle_Pl;
  }


multipole_Pk multipole_Pk_calculator::calculate_Legendre(const Mpc_units::energy& k, const Matsubara_XY& XY,
                                                         const oneloop_Pk& data, const oneloop_growth_record& gf_data,
                                                         const wiggle_Pk& Ptree)
  {
    // construct lambdas to access components of an RSD P(k) record

    auto raw_tree        = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_tree().get_raw().get_value(); };
    auto raw_13          = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_13().get_raw().get_value(); };
    auto raw_22          = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_22().get_raw().get_value(); };
    auto raw_SPT         = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_1loop_SPT().get_raw().get_value(); };
    auto raw_Z2d         = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_delta().get_raw().get_value(); };
    auto raw_Z0v         = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_Z0_v().get_raw().get_value(); };
    auto raw_Z2v         = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_v().get_raw().get_value(); };
    auto raw_Z0vd        = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_Z0_vdelta().get_raw().get_value(); };
    auto raw_Z2vd        = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_vdelta().get_raw().get_value(); };
    auto raw_Z2vv        = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_vv().get_raw().get_value(); };
    auto raw_Z2vvd       = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_vvdelta().get_raw().get_value(); };
    auto raw_Z2vvv       = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_vvv().get_raw().get_value(); };
    
    auto wiggle_tree     = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_tree().get_wiggle().get_value(); };
    auto wiggle_13       = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_13().get_wiggle().get_value(); };
    auto wiggle_22       = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_22().get_wiggle().get_value(); };
    auto wiggle_SPT      = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_1loop_SPT().get_wiggle().get_value(); };
    auto wiggle_Z2d      = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_delta().get_wiggle().get_value(); };
    auto wiggle_Z0v      = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_Z0_v().get_wiggle().get_value(); };
    auto wiggle_Z2v      = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_v().get_wiggle().get_value(); };
    auto wiggle_Z0vd     = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_Z0_vdelta().get_wiggle().get_value(); };
    auto wiggle_Z2vd     = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_vdelta().get_wiggle().get_value(); };
    auto wiggle_Z2vv     = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_vv().get_wiggle().get_value(); };
    auto wiggle_Z2vvd    = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_vvdelta().get_wiggle().get_value(); };
    auto wiggle_Z2vvv    = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_vvv().get_wiggle().get_value(); };
    
    auto nowiggle_tree   = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_tree().get_nowiggle().get_value(); };
    auto nowiggle_13     = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_13().get_nowiggle().get_value(); };
    auto nowiggle_22     = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_22().get_nowiggle().get_value(); };
    auto nowiggle_SPT    = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_1loop_SPT().get_nowiggle().get_value(); };
    auto nowiggle_Z2d    = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_delta().get_nowiggle().get_value(); };
    auto nowiggle_Z0v    = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_Z0_v().get_nowiggle().get_value(); };
    auto nowiggle_Z2v    = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_v().get_nowiggle().get_value(); };
    auto nowiggle_Z0vd   = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy3 { return data.get_Z0_vdelta().get_nowiggle().get_value(); };
    auto nowiggle_Z2vd   = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_vdelta().get_nowiggle().get_value(); };
    auto nowiggle_Z2vv   = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_vv().get_nowiggle().get_value(); };
    auto nowiggle_Z2vvd  = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_vvdelta().get_nowiggle().get_value(); };
    auto nowiggle_Z2vvv  = [&](const rsd_dd_Pk& data) -> Mpc_units::inverse_energy  { return data.get_Z2_vvv().get_nowiggle().get_value(); };
    
    // get Matsubara X+Y suppression factor (remember we have to scale up by the linear growth factor,
    // since we store just the raw integral over the early-time tree-level power spectrum)
    double Matsubara_XY = gf_data.g*gf_data.g * k*k * XY;
    
    double A_coeff = gf_data.f*(gf_data.f+2.0) * Matsubara_XY;
    double B_coeff = Matsubara_XY;
    
    multipole_Pk_calculator_impl::resum_adjuster Pk_adj(k, Matsubara_XY, gf_data, Ptree);
    multipole_Pk_calculator_impl::null_adjuster<Mpc_units::inverse_energy3> Pk_null;
    multipole_Pk_calculator_impl::null_adjuster<Mpc_units::inverse_energy> k2_Pk_null;
    
    // compute un-resummed multipole power spectra
    
    Mpc_units::inverse_energy3 P0_tree = this->decompose(raw_tree, data, multipole_Pk_calculator_impl::mu_to_ell0());
    Mpc_units::inverse_energy3 P2_tree = this->decompose(raw_tree, data, multipole_Pk_calculator_impl::mu_to_ell2());
    Mpc_units::inverse_energy3 P4_tree = this->decompose(raw_tree, data, multipole_Pk_calculator_impl::mu_to_ell4());
    
    Mpc_units::inverse_energy3 P0_13 = this->decompose(raw_13, data, multipole_Pk_calculator_impl::mu_to_ell0());
    Mpc_units::inverse_energy3 P2_13 = this->decompose(raw_13, data, multipole_Pk_calculator_impl::mu_to_ell2());
    Mpc_units::inverse_energy3 P4_13 = this->decompose(raw_13, data, multipole_Pk_calculator_impl::mu_to_ell4());
    
    Mpc_units::inverse_energy3 P0_22 = this->decompose(raw_22, data, multipole_Pk_calculator_impl::mu_to_ell0());
    Mpc_units::inverse_energy3 P2_22 = this->decompose(raw_22, data, multipole_Pk_calculator_impl::mu_to_ell2());
    Mpc_units::inverse_energy3 P4_22 = this->decompose(raw_22, data, multipole_Pk_calculator_impl::mu_to_ell4());
    
    Mpc_units::inverse_energy3 P0_SPT = this->decompose(raw_SPT, data, multipole_Pk_calculator_impl::mu_to_ell0());
    Mpc_units::inverse_energy3 P2_SPT = this->decompose(raw_SPT, data, multipole_Pk_calculator_impl::mu_to_ell2());
    Mpc_units::inverse_energy3 P4_SPT = this->decompose(raw_SPT, data, multipole_Pk_calculator_impl::mu_to_ell4());

    // compute resummed multipole power spectra
    
    Mpc_units::inverse_energy3 P0_tree_rs = this->decompose(wiggle_tree, nowiggle_tree, data, Pk_null, multipole_Pk_calculator_impl::mu_to_ell0(), multipole_Pk_calculator_impl::mu_to_ell0_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P2_tree_rs = this->decompose(wiggle_tree, nowiggle_tree, data, Pk_null, multipole_Pk_calculator_impl::mu_to_ell2(), multipole_Pk_calculator_impl::mu_to_ell2_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P4_tree_rs = this->decompose(wiggle_tree, nowiggle_tree, data, Pk_null, multipole_Pk_calculator_impl::mu_to_ell4(), multipole_Pk_calculator_impl::mu_to_ell4_expXY(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy3 P0_13_rs = this->decompose(wiggle_13, nowiggle_13, data, Pk_null, multipole_Pk_calculator_impl::mu_to_ell0(), multipole_Pk_calculator_impl::mu_to_ell0_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P2_13_rs = this->decompose(wiggle_13, nowiggle_13, data, Pk_null, multipole_Pk_calculator_impl::mu_to_ell2(), multipole_Pk_calculator_impl::mu_to_ell2_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P4_13_rs = this->decompose(wiggle_13, nowiggle_13, data, Pk_null, multipole_Pk_calculator_impl::mu_to_ell4(), multipole_Pk_calculator_impl::mu_to_ell4_expXY(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy3 P0_22_rs = this->decompose(wiggle_22, nowiggle_22, data, Pk_null, multipole_Pk_calculator_impl::mu_to_ell0(), multipole_Pk_calculator_impl::mu_to_ell0_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P2_22_rs = this->decompose(wiggle_22, nowiggle_22, data, Pk_null, multipole_Pk_calculator_impl::mu_to_ell2(), multipole_Pk_calculator_impl::mu_to_ell2_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P4_22_rs = this->decompose(wiggle_22, nowiggle_22, data, Pk_null, multipole_Pk_calculator_impl::mu_to_ell4(), multipole_Pk_calculator_impl::mu_to_ell4_expXY(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy3 P0_SPT_rs = this->decompose(wiggle_SPT, nowiggle_SPT, data, Pk_adj, multipole_Pk_calculator_impl::mu_to_ell0(), multipole_Pk_calculator_impl::mu_to_ell0_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P2_SPT_rs = this->decompose(wiggle_SPT, nowiggle_SPT, data, Pk_adj, multipole_Pk_calculator_impl::mu_to_ell2(), multipole_Pk_calculator_impl::mu_to_ell2_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P4_SPT_rs = this->decompose(wiggle_SPT, nowiggle_SPT, data, Pk_adj, multipole_Pk_calculator_impl::mu_to_ell4(), multipole_Pk_calculator_impl::mu_to_ell4_expXY(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy P0_Z2d_rs = this->decompose(wiggle_Z2d, nowiggle_Z2d, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell0(), multipole_Pk_calculator_impl::mu_to_ell0_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy P2_Z2d_rs = this->decompose(wiggle_Z2d, nowiggle_Z2d, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell2(), multipole_Pk_calculator_impl::mu_to_ell2_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy P4_Z2d_rs = this->decompose(wiggle_Z2d, nowiggle_Z2d, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell4(), multipole_Pk_calculator_impl::mu_to_ell4_expXY(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy3 P0_Z0v_rs = this->decompose(wiggle_Z0v, nowiggle_Z0v, data, Pk_null, multipole_Pk_calculator_impl::mu_to_ell0(), multipole_Pk_calculator_impl::mu_to_ell0_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P2_Z0v_rs = this->decompose(wiggle_Z0v, nowiggle_Z0v, data, Pk_null, multipole_Pk_calculator_impl::mu_to_ell2(), multipole_Pk_calculator_impl::mu_to_ell2_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P4_Z0v_rs = this->decompose(wiggle_Z0v, nowiggle_Z0v, data, Pk_null, multipole_Pk_calculator_impl::mu_to_ell4(), multipole_Pk_calculator_impl::mu_to_ell4_expXY(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy P0_Z2v_rs = this->decompose(wiggle_Z2v, nowiggle_Z2v, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell0(), multipole_Pk_calculator_impl::mu_to_ell0_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy P2_Z2v_rs = this->decompose(wiggle_Z2v, nowiggle_Z2v, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell2(), multipole_Pk_calculator_impl::mu_to_ell2_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy P4_Z2v_rs = this->decompose(wiggle_Z2v, nowiggle_Z2v, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell4(), multipole_Pk_calculator_impl::mu_to_ell4_expXY(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy3 P0_Z0vd_rs = this->decompose(wiggle_Z0vd, nowiggle_Z0vd, data, Pk_null, multipole_Pk_calculator_impl::mu_to_ell0(), multipole_Pk_calculator_impl::mu_to_ell0_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P2_Z0vd_rs = this->decompose(wiggle_Z0vd, nowiggle_Z0vd, data, Pk_null, multipole_Pk_calculator_impl::mu_to_ell2(), multipole_Pk_calculator_impl::mu_to_ell2_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy3 P4_Z0vd_rs = this->decompose(wiggle_Z0vd, nowiggle_Z0vd, data, Pk_null, multipole_Pk_calculator_impl::mu_to_ell4(), multipole_Pk_calculator_impl::mu_to_ell4_expXY(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy P0_Z2vd_rs = this->decompose(wiggle_Z2vd, nowiggle_Z2vd, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell0(), multipole_Pk_calculator_impl::mu_to_ell0_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy P2_Z2vd_rs = this->decompose(wiggle_Z2vd, nowiggle_Z2vd, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell2(), multipole_Pk_calculator_impl::mu_to_ell2_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy P4_Z2vd_rs = this->decompose(wiggle_Z2vd, nowiggle_Z2vd, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell4(), multipole_Pk_calculator_impl::mu_to_ell4_expXY(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy P0_Z2vv_rs = this->decompose(wiggle_Z2vv, nowiggle_Z2vv, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell0(), multipole_Pk_calculator_impl::mu_to_ell0_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy P2_Z2vv_rs = this->decompose(wiggle_Z2vv, nowiggle_Z2vv, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell2(), multipole_Pk_calculator_impl::mu_to_ell2_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy P4_Z2vv_rs = this->decompose(wiggle_Z2vv, nowiggle_Z2vv, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell4(), multipole_Pk_calculator_impl::mu_to_ell4_expXY(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy P0_Z2vvd_rs = this->decompose(wiggle_Z2vvd, nowiggle_Z2vvd, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell0(), multipole_Pk_calculator_impl::mu_to_ell0_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy P2_Z2vvd_rs = this->decompose(wiggle_Z2vvd, nowiggle_Z2vvd, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell2(), multipole_Pk_calculator_impl::mu_to_ell2_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy P4_Z2vvd_rs = this->decompose(wiggle_Z2vvd, nowiggle_Z2vvd, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell4(), multipole_Pk_calculator_impl::mu_to_ell4_expXY(A_coeff, B_coeff));
    
    Mpc_units::inverse_energy P0_Z2vvv_rs = this->decompose(wiggle_Z2vvv, nowiggle_Z2vvv, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell0(), multipole_Pk_calculator_impl::mu_to_ell0_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy P2_Z2vvv_rs = this->decompose(wiggle_Z2vvv, nowiggle_Z2vvv, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell2(), multipole_Pk_calculator_impl::mu_to_ell2_expXY(A_coeff, B_coeff));
    Mpc_units::inverse_energy P4_Z2vvv_rs = this->decompose(wiggle_Z2vvv, nowiggle_Z2vvv, data, k2_Pk_null, multipole_Pk_calculator_impl::mu_to_ell4(), multipole_Pk_calculator_impl::mu_to_ell4_expXY(A_coeff, B_coeff));
    
    Pk_ell P0(P0_tree, P0_tree_rs, P0_13, P0_13_rs, P0_22, P0_22_rs, P0_SPT, P0_SPT_rs,
              P0_Z2d_rs, P0_Z0v_rs, P0_Z2v_rs, P0_Z0vd_rs, P0_Z2vd_rs, P0_Z2vv_rs, P0_Z2vvd_rs, P0_Z2vvv_rs);

    Pk_ell P2(P2_tree, P2_tree_rs, P2_13, P2_13_rs, P2_22, P2_22_rs, P2_SPT, P2_SPT_rs,
              P2_Z2d_rs, P2_Z0v_rs, P2_Z2v_rs, P2_Z0vd_rs, P2_Z2vd_rs, P2_Z2vv_rs, P2_Z2vvd_rs, P2_Z2vvv_rs);

    Pk_ell P4(P4_tree, P4_tree_rs, P4_13, P4_13_rs, P4_22, P4_22_rs, P4_SPT, P4_SPT_rs,
              P4_Z2d_rs, P4_Z0v_rs, P4_Z2v_rs, P4_Z0vd_rs, P4_Z2vd_rs, P4_Z2vv_rs, P4_Z2vvd_rs, P4_Z2vvv_rs);
    
    return multipole_Pk(data.get_k_token(), data.get_Pk_token(), data.get_IR_token(), data.get_UV_token(),
                        data.get_z_token(), XY.get_IR_resum_token(), P0, P2, P4);
  }


Matsubara_XY
multipole_Pk_calculator::calculate_Matsubara_XY(const Mpc_units::energy& IR_resum, const IR_resum_token& IR_resum_tok,
                                                const wiggle_Pk& Pk_lin)
  {
    // extract database for power spectra
    const auto& raw_db = Pk_lin.get_raw_db();
    const auto& wiggle_db = Pk_lin.get_wiggle_db();
    
    // use standard clearance above lower limit of spline to avoid unwanted effects associated
    // with inaccuracies in the fit there
    const auto k_min = SPLINE_PK_DEFAULT_BOTTOM_CLEARANCE * std::min(raw_db.get_k_min(), wiggle_db.get_k_min());
    
    wiggle_Pk_nowiggle_adapter nowiggle(Pk_lin);

    // disable Cuba's built-in parallelization
    cubacores(0, multipole_Pk_calculator_impl::pcores);
    
    Mpc_units::inverse_energy2 X = this->compute_XY(IR_resum, k_min, nowiggle, multipole_Pk_calculator_impl::matsubara_X_integrand);
    Mpc_units::inverse_energy2 Y = this->compute_XY(IR_resum, k_min, nowiggle, multipole_Pk_calculator_impl::matsubara_Y_integrand);
    
    return Matsubara_XY(Pk_lin.get_token(), IR_resum_tok, X, Y);
  }


Mpc_units::inverse_energy2
multipole_Pk_calculator::compute_XY(const Mpc_units::energy& IR_resum, const Mpc_units::energy& k_min,
                                    const spline_Pk& Pk, integrand_t integrand)
  {
    cubareal integral[multipole_Pk_calculator_impl::dimensions];
    cubareal error[multipole_Pk_calculator_impl::dimensions];
    cubareal prob[multipole_Pk_calculator_impl::dimensions];
    
    int regions;
    int evaluations;
    int fail;
    
    const Mpc_units::inverse_energy qmin = 10 * Mpc_units::Mpc;
    const Mpc_units::inverse_energy qmax = 300 * Mpc_units::Mpc;
    
    std::unique_ptr<multipole_Pk_calculator_impl::integrand_data> data =
      std::make_unique<multipole_Pk_calculator_impl::integrand_data>(IR_resum, k_min, qmin, qmax, Pk);
    
    Cuhre(multipole_Pk_calculator_impl::dimensions,
          multipole_Pk_calculator_impl::components,
          integrand, data.get(),
          multipole_Pk_calculator_impl::points_per_invocation,
          this->rel_err, this->abs_err,
          multipole_Pk_calculator_impl::verbosity_none | multipole_Pk_calculator_impl::samples_last,
          multipole_Pk_calculator_impl::min_eval, multipole_Pk_calculator_impl::max_eval,
          multipole_Pk_calculator_impl::cuhre_key,
          nullptr, nullptr,
          &regions, &evaluations, &fail,
          integral, error, prob);
    
    // phase space factor is Matsubara's 1/6pi^2 multiplied by 4pi coming from the q integral angular average
    return integral[0] * Mpc_units::Mpc2 / (6.0 * M_PI * M_PI);
  }
