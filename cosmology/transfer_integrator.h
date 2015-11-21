//
// Created by David Seery on 15/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_TRANSFER_INTEGRATOR_H
#define LSSEFT_TRANSFER_INTEGRATOR_H


#include <memory>

#include "FRW_model.h"
#include "concepts/transfer_function.h"

#include "units/eV_units.h"
#include "database/tokens.h"
#include "database/z_database.h"

#include "defaults.h"

#include "boost/timer/timer.hpp"
#include "boost/serialization/serialization.hpp"


class transfer_integrator
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    transfer_integrator(double a=LSSEFT_DEFAULT_ABS_ERR, double r=LSSEFT_DEFAULT_REL_ERR);

    //! destructor is default
    ~transfer_integrator() = default;


    // TRANSFER FUNCTION

  public:

    //! integrate transfer function for a given k-mode and set of redshift samples
    transfer_function integrate(const FRW_model& model, const eV_units::energy& k, const k_token& tok,
                                std::shared_ptr<z_database>& z_db);


    // INTERNAL DATA

  private:

    // TOLERANCES

    //! required absolute error
    double abs_err;

    //! required relative error
    double rel_err;

  };


#endif //LSSEFT_TRANSFER_INTEGRATOR_H
