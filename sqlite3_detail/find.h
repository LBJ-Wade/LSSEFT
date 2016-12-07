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
#include "cosmology/concepts/Matsubara_XY.h"
#include "cosmology/concepts/power_spectrum.h"

#include "sqlite3.h"


namespace sqlite3_operations
  {
    
    //! extract one-loop growth g- and f-functions for a given set of redshifts
    std::unique_ptr<oneloop_growth>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& token,
         const z_database& z_db);
    
    //! extract filtered linear power spectrum for a given set of k-moes
    std::unique_ptr<wiggle_Pk>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const linear_Pk_token& token,
         const k_database& k_db);
    
    //! extract loop integrals for a given wavenumber, linear power spectrum, UV-cutoff and IR-cutoff combination
    std::unique_ptr<loop_integral>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
         const k_token& k, const linear_Pk_token& Pk, const IR_cutoff_token& IR_cutoff,
         const UV_cutoff_token& UV_cutoff);
    
    //! extract P(k) data for a given wavenumber, z-value, linear power spectrum, UV-cutoff and IR-cutoff combination
    std::unique_ptr<oneloop_Pk>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
             const k_token& k, const z_token& z, const linear_Pk_token& Pk_lin, const IR_cutoff_token& IR_cutoff,
             const UV_cutoff_token& UV_cutoff);
    
    //! extract Matsubara X & Y coefficient sfor a given linear power spectrum and IR resummation scale
    std::unique_ptr<Matsubara_XY>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
         const linear_Pk_token& Pk, const IR_resum_token& IR_resum);
    
  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_FIND_H
