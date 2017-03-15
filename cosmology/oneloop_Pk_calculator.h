//
// Created by David Seery on 14/11/2016.
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

#ifndef LSSEFT_ONE_LOOP_PK_CALCULATOR_H
#define LSSEFT_ONE_LOOP_PK_CALCULATOR_H


#include <list>

#include "FRW_model.h"
#include "concepts/oneloop_Pk.h"
#include "concepts/oneloop_resum_Pk.h"
#include "concepts/oneloop_growth.h"
#include "concepts/loop_integral.h"
#include "concepts/Matsubara_XY.h"
#include "concepts/power_spectrum.h"

#include "database/tokens.h"

#include "defaults.h"


class oneloop_Pk_calculator
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor is default
    oneloop_Pk_calculator() = default;
    
    //! destructor is default
    ~oneloop_Pk_calculator() = default;
    
    
    // INTERFACE
    
  public:
    
    //! calculate power spectrum
    std::list<oneloop_Pk>
    calculate_dd(const Mpc_units::energy& k, const k_token& k_tok, const oneloop_growth& gf_factors,
                     const loop_integral& loop_data, const initial_filtered_Pk& Pk_init,
                     const boost::optional<const final_filtered_Pk&>& Pk_final);
    
    //! calculate resummed power spectrum
    
    oneloop_resum_Pk
    calculate_resum_dd(const Mpc_units::energy& k, const Matsubara_XY& XY, const oneloop_Pk& data,
                       const oneloop_growth_record& gf_data, const initial_filtered_Pk& init_Pk,
                       const boost::optional<const final_filtered_Pk&>& final_Pk);

    
    // INTERNAL API
    
  private:
    
    //! compute delta-delta power spectrum and counterterm
    dd_Pk compute_dd(const Mpc_units::energy& k, const oneloop_growth_record& val, const loop_integral& loop_data,
                     const Pk_value& Ptr_init, const boost::optional<Pk_value>& Ptr_final);
    
    //! compute mu^0 coefficient of rsd delta-delta power spectrum and counterterms
    rsd_dd_Pk compute_rsd_dd_mu0(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                 const loop_integral& loop_data, const Pk_value& Ptr_init, const boost::optional<Pk_value>& Ptr_final);
    
    //! compute mu^2 coefficient of rsd delta-delta power spectrum and counterterms
    rsd_dd_Pk compute_rsd_dd_mu2(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                 const loop_integral& loop_data, const Pk_value& Ptr_init, const boost::optional<Pk_value>& Ptr_final);
    
    //! compute mu^4 coefficient of rsd delta-delta power spectrum and counterterms
    rsd_dd_Pk compute_rsd_dd_mu4(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                 const loop_integral& loop_data, const Pk_value& Ptr_init, const boost::optional<Pk_value>& Ptr_final);
    
    //! compute mu^6 coefficient of rsd delta-delta power spectrum and counterterms
    rsd_dd_Pk compute_rsd_dd_mu6(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                 const loop_integral& loop_data, const Pk_value& Ptr_init, const boost::optional<Pk_value>& Ptr_final);
    
    //! compute mu^8 coefficient of rsd delta-delta power spectrum and counterterms
    rsd_dd_Pk compute_rsd_dd_mu8(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                 const loop_integral& loop_data, const Pk_value& Ptr_init, const boost::optional<Pk_value>& Ptr_final);
    
  };


#endif //LSSEFT_ONE_LOOP_PK_CALCULATOR_H
