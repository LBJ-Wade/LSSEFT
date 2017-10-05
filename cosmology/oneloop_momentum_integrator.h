//
// Created by David Seery on 21/11/2015.
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

#ifndef LSSEFT_ONELOOP_MOMENTUM_INTEGRATOR_H
#define LSSEFT_ONELOOP_MOMENTUM_INTEGRATOR_H


#include <random>
#include <error/error_handler.h>

#include "FRW_model.h"
#include "concepts/loop_integral.h"
#include "cosmology/concepts/power_spectrum.h"

#include "database/tokens.h"

#include "defaults.h"

#include "error/error_handler.h"

#include "cuba.h"

#include "boost/timer/timer.hpp"
#include "boost/serialization/serialization.hpp"


enum class loop_integral_type { P13, P22 };


class loop_integral_params
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor
    loop_integral_params(double a_13=LSSEFT_DEFAULT_INTEGRAL_ABS_ERR_13, double r_13=LSSEFT_DEFAULT_INTEGRAL_REL_ERR_13,
                         double a_22=LSSEFT_DEFAULT_INTEGRAL_ABS_ERR_22, double r_22=LSSEFT_DEFAULT_INTEGRAL_REL_ERR_22)
      : abs_err_13(a_13),
        rel_err_13(r_13),
        abs_err_22(a_22),
        rel_err_22(r_22)
      {
      }
   
    //! destructor is default
    ~loop_integral_params() = default;
    
    
    // INTERFACE
    
  public:
    
    //! get 13 abserr
    double get_abserr_13() const { return this->abs_err_13; }
    
    //! get 13 relerr
    double get_relerr_13() const { return this->rel_err_13; }
    
    //! get 22 abserr
    double get_abserr_22() const { return this->abs_err_22; }
    
    //! get 22 relerr
    double get_relerr_22() const { return this->rel_err_22; }
    
    
    // INTERNAL DATA
  
  private:
    
    //! absolute tolerance for 13 integrals
    double abs_err_13;
    
    //! relative tolerance for 13 integrals
    double rel_err_13;
    
    //! absolute tolerance for 22 integrals
    double abs_err_22;
    
    //! relative tolerance for 22 integrals
    double rel_err_22;
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & abs_err_13;
        ar & rel_err_13;
        ar & abs_err_22;
        ar & rel_err_22;
      }
    
  };


class oneloop_momentum_integrator
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    oneloop_momentum_integrator(const loop_integral_params& p, error_handler& e);

    //! destructor is default
    ~oneloop_momentum_integrator() = default;


    // ONE-LOOP KERNEL INTEGRALS

  public:

    //! integrate one-loop kernels
    loop_integral integrate(const FRW_model& model, const loop_integral_params_token& params_tok, const Mpc_units::energy& k,
                                const k_token& k_tok, const Mpc_units::energy& UV_cutoff, const UV_cutoff_token& UV_tok,
                                const Mpc_units::energy& IR_cutoff, const IR_cutoff_token& IR_tok, const initial_filtered_Pk& Pk);

    // INTERNAL API

  private:

    //! perform a kernel integral, both raw and wiggle parts
    template <typename KernelRecord>
    bool kernel_integral(const FRW_model& model, const Mpc_units::energy& k, const Mpc_units::energy& UV_cutoff,
                         const Mpc_units::energy& IR_cutoff, const initial_filtered_Pk& Pk, integrand_t interand,
                         KernelRecord& result, loop_integral_type type, const std::string& name);
    
    //! perform a single kernel integral of either raw or wiggle type, depending on which spline is supplied
    template <typename IntegralRecord>
    bool evaluate_integral(const FRW_model& model, const Mpc_units::energy& k, const Mpc_units::energy& UV_cutoff,
                           const Mpc_units::energy& IR_cutoff, const spline_Pk& Pk, integrand_t integrand,
                           IntegralRecord& result, loop_integral_type type, const std::string& name,
                           const std::string& component);


    // INTERNAL DATA

  private:

    //! parameter block
    loop_integral_params params;
    
    //! reference to error handler agent
    error_handler& err_handler;


    // RANDOM NUMBER GENERATORS

    std::random_device random_device;
    std::mt19937       mersenne_twister;

  };


#endif //LSSEFT_ONELOOP_MOMENTUM_INTEGRATOR_H
