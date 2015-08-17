//
// Created by David Seery on 13/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_MISSING_ELEMENTS_H
#define LSSEFT_SQLITE3_MISSING_ELEMENTS_H


#include "sqlite3_policy.h"

#include "database/transaction_manager.h"
#include "database/tokens.h"
#include "database/redshift_database.h"

#include "sqlite3.h"


namespace sqlite3_operations
  {

    //! construct a database of redshifts which need to be computed for the transfer function at a given wavenumber

    //! we use a shared pointer to avoid costly duplication of the database as it gets
    //! moved around
    std::shared_ptr<redshift_database> missing_redshifts(sqlite3* db, transaction_manager& mgr,
                                                         const sqlite3_policy& policy, const FRW_model_token& model,
                                                         const wavenumber_token& k, const redshift_database& z_db,
                                                         const std::string& z_table);


    //! construct a database of redshifts which need to be computed for the one-loop growth functions

    //! we use a shared pointer to avoid costly duplication of the database as it gets
    //! moved around
    std::shared_ptr<redshift_database> missing_redshifts(sqlite3* db, transaction_manager& mgr,
                                                         const sqlite3_policy& policy, const FRW_model_token& model,
                                                         const redshift_database& z_db, const std::string& z_table);

  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_MISSING_ELEMENTS_H
