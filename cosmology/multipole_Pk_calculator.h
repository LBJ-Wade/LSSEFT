//
// Created by David Seery on 18/11/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_MULTIPOLE_PK_CALCULATOR_H
#define LSSEFT_MULTIPOLE_PK_CALCULATOR_H


#include "FRW_model.h"
#include "concepts/oneloop_growth.h"
#include "concepts/oneloop_Pk.h"
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
    multipole_Pk
    calculate_Legendre(const Mpc_units::energy& k, const Matsubara_XY& XY,
                       const oneloop_Pk& data, const oneloop_growth_record& gf_data, const wiggle_Pk& Ptree);
    
  private:
    
    //! decompose into Legendre modes using a specified decomposer functional and for a specified ell mode
    template <typename Accessor, typename Decomposer>
    auto decompose(Accessor extract, const oneloop_Pk& data, Decomposer decomp);
    
    //! decompose into Legendre modes using resummation of the wiggle part, including subtractions
    //! to prevent double-counting after resummation
    template <typename WiggleAccessor, typename NoWiggleAccessor, typename ResumAdjuster, typename RawDecomposer, typename XYDecomposer>
    auto decompose(WiggleAccessor wiggle, NoWiggleAccessor nowiggle, const oneloop_Pk& data,
                   ResumAdjuster adjust, RawDecomposer raw_decomp, XYDecomposer XY_decomp);
    
  };


#endif //LSSEFT_MULTIPOLE_PK_CALCULATOR_H
