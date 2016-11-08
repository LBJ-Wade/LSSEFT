//
// Created by David Seery on 13/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_TEMPORARY_TABLES_H
#define LSSEFT_SQLITE3_TEMPORARY_TABLES_H


#include <memory>
#include <string>

#include "sqlite3_policy.h"

#include "database/transaction_manager.h"
#include "database/tokens.h"
#include "database/z_database.h"
#include "database/k_database.h"

#include "sqlite3.h"


namespace sqlite3_operations
  {

    //! create temporary table of redshifts
    std::string z_table(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                        z_database& z_db);

    //! create temporary table of wavenumbres
    std::string k_table(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                        k_database& k_db);

    //! drop a temporary table
    void drop_temp(sqlite3* db, transaction_manager& mgr, const std::string& table);

  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_TEMPORARY_TABLES_H
