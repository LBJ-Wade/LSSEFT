//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_ONELOOP_MOMENTUM_INTEGRATOR_H
#define LSSEFT_ONELOOP_MOMENTUM_INTEGRATOR_H


#include "FRW_model.h"
#include "concepts/loop_integral.h"
#include "concepts/tree_power_spectrum.h"

#include "database/tokens.h"

#include "defaults.h"

#include "boost/timer/timer.hpp"


class oneloop_momentum_integrator
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    oneloop_momentum_integrator(double a=LSSEFT_DEFAULT_INTEGRAL_ABS_ERR, double r=LSSEFT_DEFAULT_INTEGRAL_REL_ERR);

    //! destructor is default
    ~oneloop_momentum_integrator() = default;


    // ONE-LOOP KERNEL INTEGRALS

  public:

    //! integrate one-loop kernels
    loop_integral integrate(const FRW_model& model, const eV_units::energy& k, const k_token& k_tok,
                            const eV_units::energy& UV_cutoff, const UV_token& UV_tok,
                            const eV_units::energy& IR_cutoff, const IR_token& IR_tok,
                            std::shared_ptr<tree_power_spectrum>& Pk);


    // INTERNAL API

  private:

    //! perform A integral
    eV_units::inverse_energy3 A_integral(const FRW_model& model, const eV_units::energy& k,
                                         const eV_units::energy& UV_cutoff, const eV_units::energy& IR_cutoff,
                                         std::shared_ptr<tree_power_spectrum>& Pk);


    // INTERNAL DATA

  private:

    //! absolute tolerance
    double abs_err;

    //! relative tolerance
    double rel_err;

  };


#endif //LSSEFT_ONELOOP_MOMENTUM_INTEGRATOR_H
