//
// Created by David Seery on 12/08/2015.
// --@@ // Copyright (c) 2017 University of Sussex. All rights reserved.
//
// This file is part of the Sussex Effective Field Theory for
// Large-Scale Structure platform (LSSEFT).
//
// LSSEFT is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// LSSEFT is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with LSSEFT.  If not, see <http://www.gnu.org/licenses/>.
//
// @license: GPL-2
// @contributor: David Seery <D.Seery@sussex.ac.uk>
// --@@
//

#ifndef LSSEFT_SQLITE3_WAVENUMBER_H
#define LSSEFT_SQLITE3_WAVENUMBER_H


#include <memory>

#include "database/transaction_manager.h"
#include "database/tokens.h"

#include "units/Mpc_units.h"

#include "sqlite3_policy.h"

#include "utilities.h"
#include "exceptions.h"

#include "localizations/messages.h"

#include "boost/optional.hpp"

#include "sqlite3.h"


namespace sqlite3_operations
  {

    template <typename Token>
    boost::optional<unsigned int> lookup_wavenumber(sqlite3* db, transaction_manager& mgr,
                                                    const Mpc_units::energy& k, const sqlite3_policy& policy, double tol)
      {
        assert(db != nullptr);

        double k_in_h_inv_Mpc = k * Mpc_units::Mpc;

        std::ostringstream select_stmt;
        select_stmt
        << "SELECT id FROM " << tokenization_table<Token>(policy) << " WHERE "
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


    template <typename Token>
    unsigned int insert_wavenumber(sqlite3* db, transaction_manager& mgr,
                                   const Mpc_units::energy& k, const sqlite3_policy& policy)
      {
        assert(db != nullptr);

        // get number of rows in table; this will be the identifier for the new wavenumber
        unsigned int new_id = count(db, tokenization_table<Token>(policy));

        double k_in_h_inv_Mpc = k * Mpc_units::Mpc;

        std::ostringstream insert_stmt;
        insert_stmt
        << "INSERT INTO " << tokenization_table<Token>(policy) << " VALUES (@id, @k);";

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


#endif //LSSEFT_SQLITE3_WAVENUMBER_H
