//
// Created by David Seery on 09/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#include "Matsubara_XY_calculator.h"

#include "boost/math/special_functions/bessel.hpp"


namespace Matsubara_XY_calculator_impl
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
    
    
    static int matsubara_X_integrand(const int* ndim, const cubareal* x, const int* ncomp, cubareal* f, void* userdata)
      {
        Matsubara_XY_calculator_impl::integrand_data* data = static_cast<Matsubara_XY_calculator_impl::integrand_data*>(userdata);
        
        Mpc_units::energy s = data->IR_cutoff + x[0] * data->s_range;
        Mpc_units::inverse_energy q = data->qmin + x[1] * data->q_range;
        
        Mpc_units::energy qfactor = q*q / data->q_volume;
        
        f[0] = data->jacobian * qfactor * data->Pk(s) * (1.0 - (3.0/(q*s))*boost::math::sph_bessel(1, q*s)) / Mpc_units::Mpc2;
        
        return(0);  // return value irrelevant unless = -999, which means stop integration
      }
    
    
    static int matsubara_Y_integrand(const int* ndim, const cubareal* x, const int* ncomp, cubareal* f, void* userdata)
      {
        Matsubara_XY_calculator_impl::integrand_data* data = static_cast<Matsubara_XY_calculator_impl::integrand_data*>(userdata);
        
        Mpc_units::energy s = data->IR_cutoff + x[0] * data->s_range;
        Mpc_units::inverse_energy q = data->qmin + x[1] * data->q_range;
        
        Mpc_units::energy qfactor = q*q / data->q_volume;
        
        f[0] = data->jacobian * qfactor * data->Pk(s) * boost::math::sph_bessel(2, q*s) / Mpc_units::Mpc2;
        
        return(0);  // return value irrelevant unless = -999, which means stop integration
      }
    
  }   // namespace Matsubara_XY_calculator_impl


Matsubara_XY
Matsubara_XY_calculator::calculate_Matsubara_XY(const Mpc_units::energy& IR_resum, const IR_resum_token& IR_resum_tok,
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
    cubacores(0, Matsubara_XY_calculator_impl::pcores);
    
    Mpc_units::inverse_energy2 X = this->compute_XY(IR_resum, k_min, nowiggle, Matsubara_XY_calculator_impl::matsubara_X_integrand);
    Mpc_units::inverse_energy2 Y = this->compute_XY(IR_resum, k_min, nowiggle, Matsubara_XY_calculator_impl::matsubara_Y_integrand);
    
    return Matsubara_XY(Pk_lin.get_token(), IR_resum_tok, X, Y);
  }


Mpc_units::inverse_energy2
Matsubara_XY_calculator::compute_XY(const Mpc_units::energy& IR_resum, const Mpc_units::energy& k_min,
                                    const spline_Pk& Pk, integrand_t integrand)
  {
    cubareal integral[Matsubara_XY_calculator_impl::dimensions];
    cubareal error[Matsubara_XY_calculator_impl::dimensions];
    cubareal prob[Matsubara_XY_calculator_impl::dimensions];
    
    int regions;
    int evaluations;
    int fail;
    
    const Mpc_units::inverse_energy qmin = 10 * Mpc_units::Mpc;
    const Mpc_units::inverse_energy qmax = 300 * Mpc_units::Mpc;
    
    std::unique_ptr<Matsubara_XY_calculator_impl::integrand_data> data =
      std::make_unique<Matsubara_XY_calculator_impl::integrand_data>(IR_resum, k_min, qmin, qmax, Pk);
    
    Cuhre(Matsubara_XY_calculator_impl::dimensions,
          Matsubara_XY_calculator_impl::components,
          integrand, data.get(),
          Matsubara_XY_calculator_impl::points_per_invocation,
          this->rel_err, this->abs_err,
          Matsubara_XY_calculator_impl::verbosity_none | Matsubara_XY_calculator_impl::samples_last,
          Matsubara_XY_calculator_impl::min_eval, Matsubara_XY_calculator_impl::max_eval,
          Matsubara_XY_calculator_impl::cuhre_key,
          nullptr, nullptr,
          &regions, &evaluations, &fail,
          integral, error, prob);
    
    // phase space factor is Matsubara's 1/6pi^2 multiplied by 4pi coming from the q integral angular average
    return integral[0] * Mpc_units::Mpc2 / (6.0 * M_PI * M_PI);
  }
