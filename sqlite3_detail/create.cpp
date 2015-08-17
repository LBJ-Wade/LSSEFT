//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include "utilities.h"

#include "create.h"


namespace sqlite3_operations
  {

    void create_tables(sqlite3* db, const sqlite3_policy& policy)
      {
        std::ostringstream models_stmt;
        models_stmt
          << "CREATE TABLE " << policy.FRW_model_table() << "("
          << "id INTEGER PRIMARY KEY, "
          << "omega_m DOUBLE, "
          << "omega_cc DOUBLE, "
          << "h DOUBLE, "
          << "T_CMB DOUBLE"
          << ")";

        exec(db, models_stmt.str());

        std::ostringstream z_config_stmt;
        z_config_stmt
          << "CREATE TABLE " << policy.redshift_config_table() << "("
          << "id INTEGER PRIMARY KEY, "
          << "z DOUBLE);";

        exec(db, z_config_stmt.str());

        std::ostringstream wavenumber_config_stmt;
        wavenumber_config_stmt
          << "CREATE TABLE " << policy.wavenumber_config_table() << "("
          << "id INTEGER PRIMARY KEY, "
          << "k DOUBLE);";

        exec(db, wavenumber_config_stmt.str());

        std::ostringstream transfer_stmt;
        transfer_stmt
          << "CREATE TABLE " << policy.transfer_table() << "("
          << "mid INTEGER, "
          << "kid INTEGER, "
          << "zid INTEGER, "
          << "delta_m DOUBLE, "
          << "delta_r DOUBLE, "
          << "theta_m DOUBLE, "
          << "theta_r DOUBLE, "
          << "Phi DOUBLE, "
          << "PRIMARY KEY (zid, kid), "
          << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
          << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id), "
          << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id));";

        exec(db, transfer_stmt.str());

        std::ostringstream oneloop_stmt;
        oneloop_stmt
          << "CREATE TABLE " << policy.oneloop_table() << "("
          << "mid INTEGER, "
          << "zid INTEGER, "
          << "g_linear DOUBLE, "
          << "A DOUBLE, "
          << "B DOUBLE, "
          << "D DOUBLE, "
          << "E DOUBLE, "
          << "F DOUBLE, "
          << "G DOUBLE, "
          << "oneloop DOUBLE, "
          << "PRIMARY KEY (mid, zid), "
          << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
          << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id));";

        exec(db, oneloop_stmt.str());
      }

  }   // namespace sqlite3_operations
