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

    std::shared_ptr<redshift_database> missing_redshifts(sqlite3* db, std::shared_ptr<transaction_manager>& mgr,
                                                         const sqlite3_policy& policy,
                                                         const std::shared_ptr<FRW_model_token>& model,
                                                         const std::shared_ptr<wavenumber_token>& k,
                                                         const std::shared_ptr<redshift_database>& z_db,
                                                         const std::string& z_table);

  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_MISSING_ELEMENTS_H
