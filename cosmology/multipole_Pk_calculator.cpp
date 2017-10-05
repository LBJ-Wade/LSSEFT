//
// Created by David Seery on 18/11/2016.
// --@@ // Copyright (c) 2017 University of Sussex. All rights reserved.
//
// This file is part of the Sussex Effective Field Theory for
// Large-Scale Structure platform (LSSEFT).
//
// LSSEFT is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// LSSEFT is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with LSSEFT.  If not, see <http://www.gnu.org/licenses/>.
//
// @license: GPL-2
// @contributor: David Seery <D.Seery@sussex.ac.uk>
// --@@
//

#include <cmath>

#include "multipole_Pk_calculator.h"

namespace multipole_Pk_calculator_impl
  {
    
    struct mu_to_ell0
      {
        
        constexpr double operator()(mu_power n)
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
    
    class mu_to_ell0_expXY
      {
      
      public:
        
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
        
      private:
        
        double A;
        double B;
        
      };

    struct mu_to_ell2
      {
        
        constexpr double operator()(mu_power n)
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
    
    class mu_to_ell2_expXY
      {
      
      public:
        
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
        
      private:
        
        double A;
        double B;
        
      };
    
    struct mu_to_ell4
      {
        
        constexpr double operator()(mu_power n)
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
    
    class mu_to_ell4_expXY
      {
      
      public:
        
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
        
      private:
        
        double A;
        double B;
        
      };
    
    
    // subtractions for resummation of full 1-loop power spectrum
    class resum_adjuster
      {
      
      public:

        resum_adjuster(const Mpc_units::energy& _k, double _XY, const oneloop_growth_record& _gf, const Pk_value& _Pk_tree)
          : f(_gf.f_lin),
            // the factor appearing in each subtraction is k^2 (X+Y) P_lin,w
            // here, k^2 (X+Y) is our variable XY
            // P_lin,w is the wiggle part of the density-density power spectrum
            factor(_XY * _Pk_tree.get_wiggle().get_value())
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
        
      private:
        
        double f;
        Mpc_units::inverse_energy3 factor;
        
      };
    
    
    // null subtractions
    template <typename Dimension>
    struct null_adjuster
      {

        constexpr Dimension operator()(mu_power n)
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
    
    
    // perform multipole decomposition
    template <typename MultipletType, typename PlainDecomposeMultiplet, typename ExpDecomposeMultiplet>
    class decomposer
      {
      
      public:
        
        // element_type will be Pk_element or k2_Pk_element or similar -- a container holding a value and an error
        typedef typename MultipletType::value_type::element_type element_type;
        
        // resum_type will be a container such as Pk_resum or k2_Pk_resum holding two element_types
        typedef typename MultipletType::value_type resum_type;
        
        //! constructor
        //! PlainDecomposer should be a triplet of decomposers for the Legendre modes ell = 0, ell = 2, ell = 4
        //! ExpDecomposer should be a triplet of decomposer for the Legendre modes ell = 0, ell = 2, ell = 4
        decomposer(const oneloop_Pk& d, PlainDecomposeMultiplet p, ExpDecomposeMultiplet e)
          : data(d),
            plain_decompose(std::move(p)),
            exp_decompose(std::move(e))
          {
          }
        
        // perform decomposition into a Legendre multiplet
        template <typename Accessor, typename ResumAdjuster>
        MultipletType compute(Accessor& a, ResumAdjuster& r)
          {
            // set up value groups for the raw and resummed components of each ell mode
            element_type raw_ell0 = this->decompose(a, std::get<0>(this->plain_decompose));
            element_type raw_ell2 = this->decompose(a, std::get<1>(this->plain_decompose));
            element_type raw_ell4 = this->decompose(a, std::get<2>(this->plain_decompose));
            
            element_type resum_ell0 = this->decompose(a, r, std::get<0>(this->plain_decompose), std::get<0>(this->exp_decompose));
            element_type resum_ell2 = this->decompose(a, r, std::get<1>(this->plain_decompose), std::get<1>(this->exp_decompose));
            element_type resum_ell4 = this->decompose(a, r, std::get<2>(this->plain_decompose), std::get<2>(this->exp_decompose));
            
            resum_type ell0(raw_ell0, resum_ell0);
            resum_type ell2(raw_ell2, resum_ell2);
            resum_type ell4(raw_ell4, resum_ell4);
            
            return MultipletType(ell0, ell2, ell4);
          }
        
        
        // perform decomposition without resummation
        template <typename Accessor, typename PlainDecompose>
        element_type decompose(Accessor& a, PlainDecompose& plain)
          {
            // use accessor to extract raw values
            const element_type& mu0 = a(this->data.get_dd_rsd_mu0()).get_raw();
            const element_type& mu2 = a(this->data.get_dd_rsd_mu2()).get_raw();
            const element_type& mu4 = a(this->data.get_dd_rsd_mu4()).get_raw();
            const element_type& mu6 = a(this->data.get_dd_rsd_mu6()).get_raw();
            const element_type& mu8 = a(this->data.get_dd_rsd_mu8()).get_raw();
            
            element_type raw = mu0 * plain(mu_power::mu0) + mu2 * plain(mu_power::mu2) + mu4 * plain(mu_power::mu4)
                               + mu6 * plain(mu_power::mu6) + mu8 * plain(mu_power::mu8);
            
            return raw;
          };
        
        
        // perform decomposition with resummation
        template <typename Accessor, typename ResumAdjuster, typename PlainDecompose, typename ExpDecompose>
        element_type decompose(Accessor& a, ResumAdjuster& r, PlainDecompose& plain, ExpDecompose& exp)
          {
            // use accessor to extract wiggle components
            const element_type w_mu0 = a(this->data.get_dd_rsd_mu0()).get_wiggle() + r(mu_power::mu0);
            const element_type w_mu2 = a(this->data.get_dd_rsd_mu2()).get_wiggle() + r(mu_power::mu2);
            const element_type w_mu4 = a(this->data.get_dd_rsd_mu4()).get_wiggle() + r(mu_power::mu4);
            const element_type w_mu6 = a(this->data.get_dd_rsd_mu6()).get_wiggle() + r(mu_power::mu6);
            const element_type w_mu8 = a(this->data.get_dd_rsd_mu8()).get_wiggle() + r(mu_power::mu8);
    
            // use accessor to extract no-wiggle components
            const element_type& nw_mu0 = a(this->data.get_dd_rsd_mu0()).get_nowiggle();
            const element_type& nw_mu2 = a(this->data.get_dd_rsd_mu2()).get_nowiggle();
            const element_type& nw_mu4 = a(this->data.get_dd_rsd_mu4()).get_nowiggle();
            const element_type& nw_mu6 = a(this->data.get_dd_rsd_mu6()).get_nowiggle();
            const element_type& nw_mu8 = a(this->data.get_dd_rsd_mu8()).get_nowiggle();
            
            // no-wiggle component is decomposed without resummation
            element_type nowiggle = nw_mu0 * plain(mu_power::mu0) + nw_mu2 * plain(mu_power::mu2) + nw_mu4 * plain(mu_power::mu4)
                                    + nw_mu6 * plain(mu_power::mu6) + nw_mu8 * plain(mu_power::mu8);
            
            // wiggle component is decomposed with resummation and possibly an adjustment
            element_type wiggle = w_mu0 * exp(mu_power::mu0) + w_mu2 * exp(mu_power::mu2) + w_mu4 * exp(mu_power::mu4)
                                  + w_mu6 * exp(mu_power::mu6) + w_mu8 * exp(mu_power::mu8);
            
            return nowiggle + wiggle;
          };
        
      private:
        
        const oneloop_Pk& data;
        PlainDecomposeMultiplet plain_decompose;
        ExpDecomposeMultiplet exp_decompose;
        
      };

    
    // accessors for ell0, ell2, ell4 from a Legendre_multiplet<>

    template <typename MultipletType>
    struct get_ell0
      {
        auto operator()(const MultipletType& m) -> decltype(m.get_ell0()) { return m.get_ell0(); }
      };
    
    template <typename MultipletType>
    struct get_ell2
      {
        auto operator()(const MultipletType& m) -> decltype(m.get_ell2()) { return m.get_ell2(); }
      };
    
    template <typename MultipletType>
    struct get_ell4
      {
        auto operator()(const MultipletType& m) -> decltype(m.get_ell4()) { return m.get_ell4(); }
      };
    
    
    template <typename Pk_Accessor>
    Pk_ell make_Pk_ell(const Pk_resum_multiplet& tree, const Pk_resum_multiplet& P13, const Pk_resum_multiplet& P22,
                       const Pk_resum_multiplet& PSPT, Pk_Accessor a)
      {
        return Pk_ell(a(tree), a(P13), a(P22), a(PSPT));
      };

    
  }   // namespace multipole_Pk_calculator_impl


multipole_Pk_set
multipole_Pk_calculator::calculate_Legendre(const Mpc_units::energy& k, const Matsubara_XY& XY, const oneloop_Pk_set& data,
                                            const oneloop_growth_record& Df_data, const initial_filtered_Pk& Pk_init,
                                            const boost::optional<const final_filtered_Pk&>& Pk_final)
  {
    using namespace multipole_Pk_calculator_impl;
    
    // construct lambdas to access components of an RSD P(k) record
    auto get_tree     = [](const rsd_dd_Pk& pkg) -> Pk_value     { return pkg.get_tree(); };
    auto get_13       = [](const rsd_dd_Pk& pkg) -> Pk_value     { return pkg.get_13(); };
    auto get_22       = [](const rsd_dd_Pk& pkg) -> Pk_value     { return pkg.get_22(); };
    auto get_SPT      = [](const rsd_dd_Pk& pkg) -> Pk_value     { return pkg.get_1loop_SPT(); };

    // get Matsubara X+Y suppression factor (remember we have to scale up by the square of the linear growth factor,
    // since we store just the raw integral over the early-time tree-level power spectrum)
    double Matsubara_XY = Df_data.D_lin*Df_data.D_lin * k*k * XY;
    
    double A_coeff = Df_data.f_lin*(Df_data.f_lin+2.0) * Matsubara_XY;
    double B_coeff = Matsubara_XY;
    
    // set policy objects to adjust the different mu dependences to account for resummation
    Pk_value Ptree = Df_data.D_lin*Df_data.D_lin * build_Pk_value(k, Pk_init);
    resum_adjuster Pk_adj(k, Matsubara_XY, Df_data, Ptree);
    null_adjuster<Mpc_units::inverse_energy3> Pk_null;

    // set up multiplet of plain decomposition coefficients
    mu_to_ell0 plain_ell0;
    mu_to_ell2 plain_ell2;
    mu_to_ell4 plain_ell4;
    auto plain_multiplet = std::make_tuple(plain_ell0, plain_ell2, plain_ell4);
    
    // set up multiplet of expXY decomposition coefficients
    mu_to_ell0_expXY exp_ell0(A_coeff, B_coeff);
    mu_to_ell2_expXY exp_ell2(A_coeff, B_coeff);
    mu_to_ell4_expXY exp_ell4(A_coeff, B_coeff);
    auto exp_multiplet = std::make_tuple(exp_ell0, exp_ell2, exp_ell4);

    using decompose_type = decomposer<Pk_resum_multiplet, decltype(plain_multiplet), decltype(exp_multiplet)>;

    multipole_Pk_set rval;

    auto apply = [&](std::string tag) -> void
      {
        const oneloop_Pk& Pk = data.at(tag);

        // create decomposer object
        decompose_type decomp(Pk, plain_multiplet, exp_multiplet);

        // get Legendre multiplets for each element of the power spectrum
        Pk_resum_multiplet tree = decomp.compute(get_tree, Pk_null);
        Pk_resum_multiplet P13 = decomp.compute(get_13, Pk_null);
        Pk_resum_multiplet P22 = decomp.compute(get_22, Pk_null);
        Pk_resum_multiplet PSPT = decomp.compute(get_SPT, Pk_adj);

        // reorganize these multiplets into Pk_ell containers for the ell=0, ell=2 and ell=4 modes
        Pk_ell P0 = make_Pk_ell(tree, P13, P22, PSPT, get_ell0<Pk_resum_multiplet>());
        Pk_ell P2 = make_Pk_ell(tree, P13, P22, PSPT, get_ell2<Pk_resum_multiplet>());
        Pk_ell P4 = make_Pk_ell(tree, P13, P22, PSPT, get_ell4<Pk_resum_multiplet>());

        // package everything up as as multiplet_Pk and emplace it with a matching tag
        rval.emplace(
          std::make_pair(tag,
            multipole_Pk{Pk.get_k_token(), Pk.get_growth_params(), Pk.get_loop_params(), XY.get_params_token(),
                         Pk.get_init_Pk_token(), Pk.get_final_Pk_token(), Pk.get_IR_token(), Pk.get_UV_token(),
                         Pk.get_z_token(), XY.get_IR_resum_token(), P0, P2, P4}
          )
        );
      };

#include "autogenerated/multipole_compute_stmts.cpp"

    return rval;
  }
