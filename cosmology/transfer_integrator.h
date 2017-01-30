//
// Created by David Seery on 15/08/2015.
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

#ifndef LSSEFT_TRANSFER_INTEGRATOR_H
#define LSSEFT_TRANSFER_INTEGRATOR_H


#include <memory>

#include "FRW_model.h"
#include "concepts/transfer_function.h"

#include "units/Mpc_units.h"
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
    transfer_integrator(double a= LSSEFT_DEFAULT_ODE_ABS_ERR, double r= LSSEFT_DEFAULT_ODE_REL_ERR);

    //! destructor is default
    ~transfer_integrator() = default;


    // TRANSFER FUNCTION

  public:

    //! integrate transfer function for a given k-mode and set of redshift samples
    transfer_function integrate(const FRW_model& model, const Mpc_units::energy& k, const k_token& tok,
                                const z_database& z_db);


    // INTERNAL DATA

  private:

    // TOLERANCES

    //! required absolute error
    double abs_err;

    //! required relative error
    double rel_err;

  };


#endif //LSSEFT_TRANSFER_INTEGRATOR_H
