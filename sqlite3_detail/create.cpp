//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include "utilities.h"
#include "create.h"


namespace sqlite3_operations
  {

    void create_tables(sqlite3* db)
      {
        std::ostringstream models_stmt;
        models_stmt
        << "CREATE TABLE models("
        << "id INTEGER PRIMARY KEY, "
        << "omega_m DOUBLE, "
        << "omega_cc DOUBLE, "
        << "h DOUBLE, "
        << "T_CMB DOUBLE"
        << ")";

        exec(db, models_stmt.str());
      }

  }   // namespace sqlite3_operations
