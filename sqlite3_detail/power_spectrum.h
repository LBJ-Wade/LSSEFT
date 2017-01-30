//
// Created by David Seery on 05/12/2016.
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

#ifndef LSSEFT_SQLITE3_POWER_SPECTRUM_H
#define LSSEFT_SQLITE3_POWER_SPECTRUM_H


#include "database/transaction_manager.h"
#include "database/tokens.h"

#include "cosmology/concepts/power_spectrum.h"

#include "units/Mpc_units.h"

#include "sqlite3_policy.h"

#include "utilities.h"
#include "exceptions.h"

#include "localizations/messages.h"

#include "boost/optional.hpp"

#include "sqlite3.h"


namespace sqlite3_operations
  {
    
    //! lookup ID for linear power spectrum
    boost::optional<unsigned int> lookup_Pk_linear(sqlite3* db, transaction_manager& mgr, const FRW_model_token& model,
                                                   const linear_Pk& Pk_lin, const sqlite3_policy& policy);
    
    //! insert a linear power spectrum configuraiton
    unsigned int insert_Pk_linear(sqlite3* db, transaction_manager& mgr, const FRW_model_token& model,
                                  const linear_Pk& Pk_lin, const sqlite3_policy& policy);
    
    
    
  }   // namespace sqlite3_operations

#endif //LSSEFT_SQLITE3_POWER_SPECTRUM_H
