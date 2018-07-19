//
// Created by David Seery on 09/12/2016.
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

#ifndef LSSEFT_MATSUBARA_XY_CALCULATOR_H
#define LSSEFT_MATSUBARA_XY_CALCULATOR_H


#include "cosmology/concepts/Matsubara_XY.h"
#include "cosmology/concepts/power_spectrum.h"

#include "database/tokens.h"

#include "defaults.h"

#include "cuba.h"


class MatsubaraXY_params
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! constructor
    MatsubaraXY_params(Mpc_units::inverse_energy qmn=LSSEFT_DEFAULT_RESUM_QMIN,
                       Mpc_units::inverse_energy qmx=LSSEFT_DEFAULT_RESUM_QMAX,
                       double a=LSSEFT_DEFAULT_INTEGRAL_ABS_ERR_22, double r=LSSEFT_DEFAULT_INTEGRAL_REL_ERR_22)
      : abs_err(a),
        rel_err(r),
        qmin(qmn),
        qmax(qmx)
      {
      }
    
    //! destructor is default
    ~MatsubaraXY_params() = default;
    
    
    // INTERFACE
  
  public:
    
    //! get abserr
    double get_abserr() const { return this->abs_err; }
    
    //! get relerr
    double get_relerr() const { return this->rel_err; }
    
    //! get qmin
    const Mpc_units::inverse_energy& get_qmin() const { return this->qmin; }
    
    //! get qmax
    const Mpc_units::inverse_energy& get_qmax() const { return this->qmax; }
    
    
    // INTERNAL DATA
  
  private:
    
    //! absolute tolerance
    double abs_err;
    
    //! relative tolerance
    double rel_err;
    
    //! minimum scale to use in averaging
    Mpc_units::inverse_energy qmin;
    
    //! maximum scale to use in averaging
    Mpc_units::inverse_energy qmax;
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & abs_err;
        ar & rel_err;
        ar & qmin;
        ar & qmax;
      }
    
  };


class Matsubara_XY_calculator
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor
    Matsubara_XY_calculator(const MatsubaraXY_params& p)
      : params(p)
      {
      }
    
    //! destructor is default
    ~Matsubara_XY_calculator() = default;
    
    
    // INTERFACE
    
  public:
    
    //! calculate Matsubara X & Y coefficients
    Matsubara_XY
    calculate_Matsubara_XY(const Mpc_units::energy& IR_resum, const IR_resum_token& IR_resum_tok,
                           const initial_filtered_Pk& Pk_lin,
                           const MatsubaraXY_params_token& params_tok);

    
    // INTERNAL API
    
  private:

    //! compute integrals for Matsubara X & Y factors
    Mpc_units::inverse_energy2
    compute_XY(const Mpc_units::energy& IR_resum, const Mpc_units::energy& k_min,
               const generic_Pk& Pk, integrand_t integrand);
    
    
    // INTERNAL DATA
  
  private:
    
    //! parameter block
    MatsubaraXY_params params;
    
  };


#endif //LSSEFT_MATSUBARA_XY_CALCULATOR_H
