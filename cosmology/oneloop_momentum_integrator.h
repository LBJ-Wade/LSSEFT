//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_ONELOOP_MOMENTUM_INTEGRATOR_H
#define LSSEFT_ONELOOP_MOMENTUM_INTEGRATOR_H


#include <random>

#include "FRW_model.h"
#include "concepts/loop_integral.h"
#include "cosmology/concepts/power_spectrum.h"

#include "database/tokens.h"

#include "defaults.h"

#include "cuba.h"

#include "boost/timer/timer.hpp"


enum class loop_integral_type { P13, P22 };


class oneloop_momentum_integrator
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    oneloop_momentum_integrator(double a_13=LSSEFT_DEFAULT_INTEGRAL_ABS_ERR_13, double r_13=LSSEFT_DEFAULT_INTEGRAL_REL_ERR_13,
                                double a_22=LSSEFT_DEFAULT_INTEGRAL_ABS_ERR_22, double r_22=LSSEFT_DEFAULT_INTEGRAL_REL_ERR_22);

    //! destructor is default
    ~oneloop_momentum_integrator() = default;


    // ONE-LOOP KERNEL INTEGRALS

  public:

    //! integrate one-loop kernels
    loop_integral integrate(const FRW_model& model, const Mpc_units::energy& k, const k_token& k_tok,
                            const Mpc_units::energy& UV_cutoff, const UV_cutoff_token& UV_tok,
                            const Mpc_units::energy& IR_cutoff, const IR_cutoff_token& IR_tok,
                            const wiggle_Pk& Pk);

    //! output integrands for inspection
    void write_integrands(const FRW_model& model, const Mpc_units::energy& k,
                          const Mpc_units::energy& UV_cutoff, const Mpc_units::energy& IR_cutoff,
                          const spline_Pk& Pk, unsigned int Npoints);

    // INTERNAL API

  private:

    //! perform a kernel integral, both raw and wiggle parts
    template <typename KernelRecord>
    bool kernel_integral(const FRW_model& model, const Mpc_units::energy& k, const Mpc_units::energy& UV_cutoff,
                         const Mpc_units::energy& IR_cutoff, const wiggle_Pk& Pk, integrand_t interand,
                         KernelRecord& result, loop_integral_type type, const std::string& name);
    
    //! perform a single kernel integral of either raw or wiggle type, depending on which spline is supplied
    template <typename IntegralRecord>
    bool evaluate_integral(const FRW_model& model, const Mpc_units::energy& k, const Mpc_units::energy& UV_cutoff,
                           const Mpc_units::energy& IR_cutoff, const spline_Pk& Pk, integrand_t integrand,
                           IntegralRecord& result, loop_integral_type type, const std::string& name,
                           const std::string& component);


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


    // RANDOM NUMBER GENERATORS

    std::random_device random_device;
    std::mt19937       mersenne_twister;

  };


#endif //LSSEFT_ONELOOP_MOMENTUM_INTEGRATOR_H
