//
// Created by David Seery on 17/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_ONELOOP_INTEGRATOR_H
#define LSSEFT_ONELOOP_INTEGRATOR_H


#include <memory>

#include "FRW_model.h"
#include "concepts/oneloop_growth.h"

#include "database/tokens.h"
#include "database/z_database.h"

#include "defaults.h"

#include "boost/timer/timer.hpp"


struct growth_integrator_data
  {
    
    growth_integrator_data(std::unique_ptr<oneloop_growth> c, boost::timer::nanosecond_type t, size_t s)
      : container(std::move(c)),
        time(t),
        steps(s)
      {
      }
    
    std::unique_ptr<oneloop_growth> container;
    boost::timer::nanosecond_type time;
    size_t steps;
    
  };


class oneloop_growth_integrator
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    oneloop_growth_integrator(double a= LSSEFT_DEFAULT_ODE_ABS_ERR, double r= LSSEFT_DEFAULT_ODE_REL_ERR);

    //! destructor is default
    ~oneloop_growth_integrator() = default;


    // ONE-LOOP GROWTH FACTORS

  public:

    //! integrate one-loop growth factors for a given set of redshift samples
    growth_integrator_data integrate(const FRW_model& model, z_database& z_db);


    // INTERNAL DATA

  private:

    // TOLERANCES

    //! required absolute error
    double abs_err;

    //! required relative error
    double rel_err;

  };


#endif //LSSEFT_ONELOOP_INTEGRATOR_H
