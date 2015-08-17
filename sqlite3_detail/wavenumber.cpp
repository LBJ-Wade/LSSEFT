//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include <sstream>
#include <assert.h>

#include "utilities.h"
#include "wavenumber.h"

#include "exceptions.h"
#include "localizations/messages.h"


namespace sqlite3_operations
  {

    boost::optional<unsigned int> lookup_wavenumber(sqlite3* db, transaction_manager& mgr,
                                                    const eV_units::energy& k, const sqlite3_policy& policy,
                                                    double tol)
      {
        assert(db != nullptr);

        double k_in_h_inv_Mpc = k * eV_units::Mpc;

        std::ostringstream select_stmt;
        select_stmt
          << "SELECT id FROM " << policy.wavenumber_config_table() << " WHERE "
          << "ABS((k-@k)/k)<@tol;";

        // prepare SQL statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, select_stmt.str().c_str(), select_stmt.str().length()+1, &stmt, nullptr));

        // bind values to the parameters
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@tol"), tol));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@k"), k_in_h_inv_Mpc));

        // execute statement and step through results
        int status = 0;
        boost::optional<unsigned int> id = boost::none;
        while((status = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(status == SQLITE_ROW)
              {
                if(id) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_MULTIPLE_WAVENUMBERS);
                id = static_cast<unsigned int>(sqlite3_column_int(stmt, 0));
              }
          }

        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));

        return(id);
      }


    unsigned int insert_wavenumber(sqlite3* db, transaction_manager& mgr,
                                   const eV_units::energy& k, const sqlite3_policy& policy)
      {
        assert(db != nullptr);

        // get number of rows in table; this will be the identifier for the new wavenumber
        unsigned int new_id = count(db, policy.wavenumber_config_table());

        double k_in_h_inv_Mpc = k * eV_units::Mpc;

        std::ostringstream insert_stmt;
        insert_stmt
          << "INSERT INTO " << policy.wavenumber_config_table() << " VALUES (@id, @k);";

        // prepare SQL statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));

        // bind values to the parameters
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@id"), new_id));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@k"), k_in_h_inv_Mpc));

        // perform insertion
        check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_WAVENUMBER_FAIL, SQLITE_DONE);

        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));

        return(new_id);
      }


  }   // namespace sqlite3_operations