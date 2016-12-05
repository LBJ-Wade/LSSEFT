//
// Created by David Seery on 05/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
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
                                                   const linear_power_spectrum& Pk_lin, const sqlite3_policy& policy);
    
    //! insert a linear power spectrum configuraiton
    unsigned int insert_Pk_linear(sqlite3* db, transaction_manager& mgr, const FRW_model_token& model,
                                  const linear_power_spectrum& Pk_lin, const sqlite3_policy& policy);
    
    
    
  }   // namespace sqlite3_operations

#endif //LSSEFT_SQLITE3_POWER_SPECTRUM_H
