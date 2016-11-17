//
// Created by David Seery on 10/11/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_FIND_H
#define LSSEFT_SQLITE3_FIND_H


#include "sqlite3_policy.h"

#include "database/transaction_manager.h"
#include "database/tokens.h"
#include "database/z_database.h"
#include "database/k_database.h"

#include "cosmology/concepts/oneloop_growth.h"
#include "cosmology/concepts/loop_integral.h"
#include "cosmology/concepts/oneloop_Pk.h"

#include "sqlite3.h"


namespace sqlite3_operations
  {
    
    //! extract one-loop growth g- and f-functions for a given set of redshifts
    oneloop_growth find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                        const FRW_model_token& token, z_database& z_db);
    
    //! extract loop integrals for a given wavenumber, UV-cutoff and IR-cutoff combination
    loop_integral find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                       const FRW_model_token& model, const k_token& k,
                       const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff);
    
    //! extract P(k) data for a given wavenumber, z-value, UV-cutoff and IR-cutoff combination
    oneloop_Pk find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                    const FRW_model_token& model, const k_token& k, const z_token& z,
                    const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff);
    
  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_FIND_H
