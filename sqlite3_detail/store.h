//
// Created by David Seery on 17/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_STORE_H
#define LSSEFT_SQLITE3_STORE_H


#include "database/transaction_manager.h"

#include "cosmology/concepts/transfer_function.h"
#include "cosmology/concepts/oneloop.h"

#include "sqlite3_policy.h"

#include "sqlite3.h"


namespace sqlite3_operations
  {

    //! store a transfer function sample
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const transfer_function& sample);

    //! store a one-loop growth factor sample
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const oneloop& sample);

  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_STORE_H
