//
// Created by David Seery on 13/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include <iostream>
#include <sstream>

#include "utilities.h"
#include "temporary_tables.h"

#include "sqlite3_defaults.h"

#include "units/Mpc_units.h"

#include "localizations/messages.h"


namespace sqlite3_operations
  {

    static unsigned int table_number = 0;

    std::string table_name(const sqlite3_policy& policy)
      {
        std::ostringstream name;

        name << policy.temp_table() << "_" << table_number++;
        return(name.str());
      }

    std::string z_table(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                        z_database& z_db)
      {
        assert(db != nullptr);

        // get new temporary table name
        std::string name = table_name(policy);

        // create table
        std::ostringstream create_stmt;
        create_stmt
          << "CREATE TEMPORARY TABLE " << name << "("
          << "id INTEGER, "
          << "z DOUBLE);";

        exec(db, create_stmt.str());

        // set up SQL insert statement
        std::ostringstream insert_stmt;
        insert_stmt
          << "INSERT INTO temp." << name << " VALUES (@id, @z);";

        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));

        // loop through records in the database, writing entries to the table
        for(z_database::record_iterator t = z_db.record_begin(); t != z_db.record_end(); ++t)
          {
            // bind parameter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@id"), t->get_token().get_id()));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@z"), *(*t)));

            // write this row
            check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_TEMPORARY_REDSHIFT, SQLITE_DONE);

            // release bindings and reset statement for next row
            check_stmt(db, sqlite3_clear_bindings(stmt));
            check_stmt(db, sqlite3_reset(stmt));
          }

        // finalize statement before exiting
        check_stmt(db, sqlite3_finalize(stmt));

        return(name);
      }


    std::string k_table(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                        k_database& k_db)
      {
        assert(db != nullptr);

        // get new temporary table name
        std::string name = table_name(policy);

        // create table
        std::ostringstream create_stmt;
        create_stmt
          << "CREATE TEMPORARY TABLE " << name << "("
          << "id INTEGER, "
          << "k DOUBLE);";

        exec(db, create_stmt.str());

        // set up SQL insert statement
        std::ostringstream insert_stmt;
        insert_stmt
          << "INSERT INTO temp." << name << " VALUES (@id, @k);";

        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));

        // loop through records in the database, writing entries to the table
        for(k_database::record_iterator t = k_db.record_begin(); t != k_db.record_end(); ++t)
          {
            double k_in_h_inv_Mpc = *(*t) * Mpc_units::Mpc;

            // bind parameter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@id"), t->get_token().get_id()));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@k"), k_in_h_inv_Mpc));

            // write this row
            check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_TEMPORARY_WAVENUMBER, SQLITE_DONE);

            // release bindings and reset statement for next row
            check_stmt(db, sqlite3_clear_bindings(stmt));
            check_stmt(db, sqlite3_reset(stmt));
          }

        // finalize statement before exiting
        check_stmt(db, sqlite3_finalize(stmt));

        return(name);
      }


    void drop_temp(sqlite3* db, transaction_manager& mgr, const std::string& table)
      {
        assert(db != nullptr);

        std::ostringstream drop_stmt;
        drop_stmt << "DROP TABLE temp." << table;

        exec(db, drop_stmt.str());
      }

  }   // namespace sqlite3_operations
