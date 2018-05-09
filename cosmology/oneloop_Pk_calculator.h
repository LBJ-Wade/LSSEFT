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
#include <map>

#include "FRW_model.h"
#include "cosmology/concepts/oneloop_Pk.h"
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
    std::list<oneloop_Pk_set>
    calculate_Pk(const Mpc_units::energy& k, const k_token& k_tok, const oneloop_growth& gf_factors,
                 const loop_integral& loop_data, const initial_filtered_Pk& Pk_init,
                 const boost::optional<const final_filtered_Pk&>& Pk_final);
    
    //! calculate resummed power spectrum

  };


#endif //LSSEFT_ONE_LOOP_PK_CALCULATOR_H
