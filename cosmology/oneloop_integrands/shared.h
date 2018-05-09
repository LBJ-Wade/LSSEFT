//
// Created by David Seery on 15/11/2016.
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

#ifndef LSSEFT_SHARED_H
#define LSSEFT_SHARED_H


#include <cmath>

#include "cosmology/FRW_model.h"
#include "cosmology/concepts/power_spectrum.h"

#include "units/Mpc_units.h"


namespace oneloop_momentum_impl
  {
    
    constexpr unsigned int dimensions            = 2;   // integrals are dq dx, even where the x integral is trivial; Cuhre seems to have issues in one dimension
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
        
        integrand_data(const FRW_model& m, const Mpc_units::energy& _k, const Mpc_units::energy& UV, const Mpc_units::energy& IR,
                       const spline_Pk& _Pk)
          : model(m),
            k(_k),
            UV_cutoff(UV),
            IR_cutoff(IR),
            Pk(_Pk),
            // Jacobian for a dq dx type integral
            // from dx there is a 2, and from dq we get UV-IR as usual
            // the difference to d^3 q is that there is no integral over theta, so no 2pi to account for
            jacobian_dqdx(2.0 * (UV_cutoff - IR_cutoff)),
            // Jacobian for a dq type integral (ie., over q alone) -- no angular directions
            jacobian_dq(UV_cutoff - IR_cutoff),
            q_range(UV_cutoff - IR_cutoff),
            k_sq(k*k)
          {
          }
        
        const FRW_model& model;
        const Mpc_units::energy& k;
        const Mpc_units::energy& UV_cutoff;
        const Mpc_units::energy& IR_cutoff;
        const spline_Pk& Pk;
        
        Mpc_units::energy  jacobian_dqdx;
        Mpc_units::energy  jacobian_dq;
        Mpc_units::energy  q_range;
        Mpc_units::energy2 k_sq;
      };
    
  }   // namespace oneloop_momentum_impl

#endif //LSSEFT_SHARED_H
