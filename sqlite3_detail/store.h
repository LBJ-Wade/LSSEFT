//
// Created by David Seery on 17/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_STORE_H
#define LSSEFT_SQLITE3_STORE_H


#include "database/transaction_manager.h"

#include "cosmology/concepts/transfer_function.h"
#include "cosmology/concepts/filtered_Pk.h"
#include "cosmology/concepts/oneloop_growth.h"
#include "cosmology/concepts/loop_integral.h"
#include "cosmology/concepts/oneloop_Pk.h"
#include "cosmology/concepts/oneloop_resum_Pk.h"
#include "cosmology/concepts/multipole_Pk.h"
#include "cosmology/concepts/Matsubara_XY.h"

#include "sqlite3_policy.h"

#include "sqlite3.h"


namespace sqlite3_operations
  {

    //! store a transfer function sample
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const transfer_function& sample);
    
    //! store a filtered power spectrum sample
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token&, const filtered_Pk& sample);

    //! store a one-loop growth factor sample
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const oneloop_growth& sample);

    //! store a loop momentum integral sample
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const loop_integral& sample);
    
    //! store Matsubara X & Y coefficients
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const Matsubara_XY& sample);
    
    //! store a one-loop Pk sample
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const oneloop_Pk& sample);
    
    //! store a resummed one-loop Pk sample
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const oneloop_resum_Pk& sample);
    
    //! store a multipole Pk sample
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const multipole_Pk& sample);

  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_STORE_H
