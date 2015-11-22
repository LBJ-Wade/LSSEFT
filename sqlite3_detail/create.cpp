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

        std::ostringstream IR_config_stmt;
        IR_config_stmt
          << "CREATE TABLE " << policy.IR_config_table() << "("
          << "id INTEGER PRIMARY KEY, "
          << "k DOUBLE);";

        exec(db, IR_config_stmt.str());

        std::ostringstream UV_config_stmt;
        UV_config_stmt
          << "CREATE TABLE " << policy.UV_config_table() << "("
          << "id INTEGER PRIMARY KEY, "
          << "k DOUBLE);";

        exec(db, UV_config_stmt.str());

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

        std::ostringstream oneloop_growth_stmt;
        oneloop_growth_stmt
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
          << "PRIMARY KEY (mid, zid), "
          << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
          << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id));";

        exec(db, oneloop_growth_stmt.str());

        std::ostringstream oneloop_momentum_stmt;
        oneloop_momentum_stmt
          << "CREATE TABLE " << policy.loop_momentum_table() << "("
          << "mid INTEGER, "
          << "kid INTEGER, "
          << "IR_id INTEGER, "
          << "UV_id INTEGER, "
          << "AA DOUBLE, "
          << "AB DOUBLE, "
          << "BB DOUBLE, "
          << "D DOUBLE, "
          << "E DOUBLE, "
          << "F DOUBLE, "
          << "G DOUBLE, "
          << "AA_regions INTEGER, "
          << "AB_regions INTEGER, "
          << "BB_regions INTEGER, "
          << "D_regions INTEGER, "
          << "E_regions INTEGER, "
          << "F_regions INTEGER, "
          << "G_regions INTEGER, "
          << "AA_evals INTEGER, "
          << "AB_evals INTEGER, "
          << "BB_evals INTEGER, "
          << "D_evals INTEGER, "
          << "E_evals INTEGER, "
          << "F_evals INTEGER, "
          << "G_evals INTEGER, "
          << "AA_err DOUBLE, "
          << "AB_err DOUBLE, "
          << "BB_err DOUBLE, "
          << "D_err DOUBLE, "
          << "E_err DOUBLE, "
          << "F_err DOUBLE, "
          << "G_err DOUBLE, "
          << "integration_time INTEGER, "
          << "PRIMARY KEY (mid, kid, IR_id, UV_id), "
          << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
          << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id), "
          << "FOREIGN KEY (IR_id) REFERENCES " << policy.IR_config_table() << "(id), "
          << "FOREIGN KEY (UV_id) REFERENCES " << policy.UV_config_table() << "(id));";

        exec(db, oneloop_momentum_stmt.str());
      }

  }   // namespace sqlite3_operations
