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

        std::ostringstream delta_m_stmt;
        delta_m_stmt
          << "CREATE TABLE " << policy.delta_m_table() << "("
          << "zid INTEGER, "
          << "kid INTEGER, "
          << "value DOUBLE, "
          << "PRIMARY KEY (zid, kid), "
          << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id), "
          << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id));";

        exec(db, delta_m_stmt.str());

        std::ostringstream delta_r_stmt;
        delta_r_stmt
        << "CREATE TABLE " << policy.delta_r_table() << "("
        << "zid INTEGER, "
        << "kid INTEGER, "
        << "value DOUBLE, "
        << "PRIMARY KEY (zid, kid), "
        << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id), "
        << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id));";

        exec(db, delta_r_stmt.str());

        std::ostringstream theta_m_stmt;
        theta_m_stmt
        << "CREATE TABLE " << policy.theta_m_table() << "("
        << "zid INTEGER, "
        << "kid INTEGER, "
        << "value DOUBLE, "
        << "PRIMARY KEY (zid, kid), "
        << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id), "
        << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id));";

        exec(db, theta_m_stmt.str());

        std::ostringstream theta_r_stmt;
        theta_r_stmt
        << "CREATE TABLE " << policy.theta_r_table() << "("
        << "zid INTEGER, "
        << "kid INTEGER, "
        << "value DOUBLE, "
        << "PRIMARY KEY (zid, kid), "
        << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id), "
        << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id));";

        exec(db, theta_r_stmt.str());

        std::ostringstream Phi_stmt;
        Phi_stmt
        << "CREATE TABLE " << policy.Phi_table() << "("
        << "zid INTEGER, "
        << "kid INTEGER, "
        << "value DOUBLE, "
        << "PRIMARY KEY (zid, kid), "
        << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id), "
        << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id));";

        exec(db, Phi_stmt.str());

        std::ostringstream oneloop_stmt;
        oneloop_stmt
        << "CREATE TABLE " << policy.oneloop_table() << "("
        << "zid INTEGER PRIMARY KEY, "
        << "A DOUBLE, "
        << "B DOUBLE, "
        << "D DOUBLE, "
        << "E DOUBLE, "
        << "F DOUBLE, "
        << "G DOUBLE, "
        << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id));";

        exec(db, oneloop_stmt.str());
      }

  }   // namespace sqlite3_operations
