//
// Created by David Seery on 13/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_MISSING_ELEMENTS_H
#define LSSEFT_SQLITE3_MISSING_ELEMENTS_H


#include "sqlite3_policy.h"

#include "database/transaction_manager.h"
#include "database/tokens.h"
#include "database/z_database.h"
#include "database/k_database.h"
#include "database/data_manager_impl.h"

#include "sqlite3.h"


namespace sqlite3_operations
  {

    //! construct a database of redshifts which need to be computed for the transfer function at a given wavenumber
    //! ownership of the resulting database is transferred via std::unique_ptr<>
    std::unique_ptr<z_database> missing_transfer_redshifts(sqlite3* db, transaction_manager& mgr,
                                                           const sqlite3_policy& policy, const FRW_model_token& model,
                                                           const k_token& k, const z_database& z_db,
                                                           const std::string& z_table);


    //! construct a database of redshifts which need to be computed for the one-loop growth functions
    //! ownership of the resulting database is transferred via std::unique_ptr<>
    std::unique_ptr<z_database> missing_oneloop_redshifts(sqlite3* db, transaction_manager& mgr,
                                                          const sqlite3_policy& policy, const FRW_model_token& model,
                                                          const z_database& z_db, const std::string& z_table);


    //! process a list of configurations for loop momentum integrals;
    //! we detect which ones are already present in the database and mark them not to be computed
    void missing_loop_momentum(sqlite3* db, transaction_manager& mgr,
                               const sqlite3_policy& policy, const FRW_model_token& model, const k_database& k_db,
                               std::list<data_manager_impl::loop_momentum_configuration>& combinations);

  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_MISSING_ELEMENTS_H
