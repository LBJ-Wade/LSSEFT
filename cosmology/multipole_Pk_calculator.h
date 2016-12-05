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
#include "concepts/Matsubara_A.h"
#include "cosmology/concepts/power_spectrum.h"

#include "database/tokens.h"

#include "defaults.h"


enum class mu_power { mu0, mu2, mu4, mu6, mu8 };


class multipole_Pk_calculator
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor
    multipole_Pk_calculator(double r = LSSEFT_DEFAULT_INTEGRAL_REL_ERR_22,
                            double a = LSSEFT_DEFAULT_INTEGRAL_ABS_ERR_22)
      : rel_err(r),
        abs_err(a)
      {
      }
    
    //! destructor is default
    ~multipole_Pk_calculator() = default;
    
    
    // INTERFACE
    
  public:
    
    //! calculate power spectra decomposition into Legendre modes
    multipole_Pk
    calculate_Legendre(const Mpc_units::energy& k, const Matsubara_A& A,
                       const oneloop_Pk& data, const oneloop_growth_record& gf_data, const tree_Pk& Ptree);
    
    //! calculate Matsubara-A coefficient
    Matsubara_A
    calculate_Matsubara_A(const Mpc_units::energy& IR_resum, const IR_resum_token& IR_resum_tok,
                          const tree_Pk& Ptree);
    
  private:
    
    //! decompose into Legendre modes using a specified decomposer functional and for a specified ell mode
    template <typename Accessor, typename Decomposer>
    decltype(std::declval<Accessor>()(std::declval<const rsd_dd_Pk&>()))
    decompose(Accessor extract, const oneloop_Pk& data, Decomposer decomp);
    
    //! decompose into Legendre modes using a specified decomposer, including the additive correction
    //! to prevent double-counting after resummation
    template <typename Decomposer>
    Mpc_units::inverse_energy3
    decompose_1loop_resummed(const oneloop_Pk& data, double Matsubara_A, const oneloop_growth_record& gf_data,
                             Decomposer decomp, const Mpc_units::energy& k, const tree_Pk& Ptree);
    
    // INTERNAL DATA
    
  private:
    
    //! relative error
    double rel_err;
    
    //! absolute error
    double abs_err;
    
  };


#endif //LSSEFT_MULTIPOLE_PK_CALCULATOR_H
