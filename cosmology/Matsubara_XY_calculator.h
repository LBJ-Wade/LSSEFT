//
// Created by David Seery on 09/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_MATSUBARA_XY_CALCULATOR_H
#define LSSEFT_MATSUBARA_XY_CALCULATOR_H


#include "cosmology/concepts/Matsubara_XY.h"
#include "cosmology/concepts/power_spectrum.h"

#include "database/tokens.h"

#include "defaults.h"

#include "cuba.h"

class Matsubara_XY_calculator
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor
    Matsubara_XY_calculator(double r = LSSEFT_DEFAULT_INTEGRAL_REL_ERR_22,
                            double a = LSSEFT_DEFAULT_INTEGRAL_ABS_ERR_22)
      : rel_err(std::abs(r)),
        abs_err(std::abs(a))
      {
      }
    
    //! destructor is default
    ~Matsubara_XY_calculator() = default;
    
    
    // INTERFACE
    
  public:
    
    //! calculate Matsubara X & Y coefficients
    Matsubara_XY
    calculate_Matsubara_XY(const Mpc_units::energy& IR_resum, const IR_resum_token& IR_resum_tok,
                           const wiggle_Pk& Pk_lin);

    
    // INTERNAL API
    
  private:

    //! compute integrals for Matsubara X & Y factors
    Mpc_units::inverse_energy2
    compute_XY(const Mpc_units::energy& IR_resum, const Mpc_units::energy& k_min,
               const spline_Pk& Pk, integrand_t integrand);
    
    
    // INTERNAL DATA
  
  private:
    
    //! relative error
    double rel_err;
    
    //! absolute error
    double abs_err;
    
  };


#endif //LSSEFT_MATSUBARA_XY_CALCULATOR_H
