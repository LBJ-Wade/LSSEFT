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

#ifndef LSSEFT_MULTIPOLE_PK_CALCULATOR_H
#define LSSEFT_MULTIPOLE_PK_CALCULATOR_H


#include "FRW_model.h"
#include "concepts/oneloop_growth.h"
#include "cosmology/concepts/oneloop_Pk.h"
#include "concepts/multipole_Pk.h"
#include "cosmology/concepts/Matsubara_XY.h"
#include "cosmology/concepts/power_spectrum.h"

#include "database/tokens.h"

#include "defaults.h"


enum class mu_power { mu0, mu2, mu4, mu6, mu8 };


class multipole_Pk_calculator
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor
    multipole_Pk_calculator() = default;
    
    //! destructor is default
    ~multipole_Pk_calculator() = default;
    
    
    // INTERFACE
    
  public:
    
    //! calculate power spectra decomposition into Legendre modes
    multipole_Pk_set
    calculate_Legendre(const Mpc_units::energy& k, const Matsubara_XY& XY, const oneloop_Pk_set& data,
                       const oneloop_growth_record& Df_data, const initial_filtered_Pk& Pk_init,
                       const boost::optional<const final_filtered_Pk&>& Pk_final);

    //! calculate counterterm decomposition into Legendre modes
    multipole_counterterm_set
    calculate_counterterms(const Mpc_units::energy& k, const k_token& k_token, const IR_cutoff_token& IR_tok,
                               const UV_cutoff_token& UV_tok, const z_token& z_tok, const growth_params_token& g_tok,
                               const Matsubara_XY& XY, const oneloop_growth_record& Df_data, const initial_filtered_Pk& Pk_init,
                               const boost::optional<const final_filtered_Pk&>& Pk_final);
    
  private:
    
    //! decompose into Legendre modes using a specified decomposer functional and for a specified ell mode
    template <typename Accessor, typename Decomposer>
    auto decompose(Accessor extract, const oneloop_Pk_set& data, Decomposer decomp);
    
    //! decompose into Legendre modes using resummation of the wiggle part, including subtractions
    //! to prevent double-counting after resummation
    template <typename WiggleAccessor, typename NoWiggleAccessor, typename ResumAdjuster, typename RawDecomposer, typename XYDecomposer>
    auto decompose(WiggleAccessor wiggle, NoWiggleAccessor nowiggle, const oneloop_Pk_set& data,
                   ResumAdjuster adjust, RawDecomposer raw_decomp, XYDecomposer XY_decomp);
    
  };


#endif //LSSEFT_MULTIPOLE_PK_CALCULATOR_H
