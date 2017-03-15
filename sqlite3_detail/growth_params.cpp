//
// Created by David Seery on 14/03/2017.
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

#include "growth_params.h"

#include "database/tokens.h"

#include "utilities.h"

#include "exceptions.h"
#include "localizations/messages.h"


namespace sqlite3_operations
  {
    
    boost::optional<unsigned int>
    lookup_growth_params(sqlite3* db, transaction_manager& mgr, const growth_params& data,
                         const sqlite3_policy& policy, double tol)
      {
        assert(db != nullptr);
        
        std::ostringstream select_stmt;
        select_stmt
          << "SELECT id FROM " << tokenization_table<growth_params_token>(policy) << " WHERE "
          << "ABS((abserr-@abs)/abserr)<@tol "
          << "AND ABS((relerr-@rel)/relerr)<@tol;";
        
        // prepare SQL statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, select_stmt.str().c_str(), select_stmt.str().length()+1, &stmt, nullptr));
        
        // bind values to the parameters
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@tol"), tol));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@abs"), data.get_abserr()));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@rel"), data.get_relerr()));
        
        // execute statement and step through results
        int status = 0;
        boost::optional<unsigned int> id = boost::none;
        while((status = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(status == SQLITE_ROW)
              {
                if(id) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_MULTIPLE_GROWTH_PARAMS);
                id = static_cast<unsigned int>(sqlite3_column_int(stmt, 0));
              }
          }
        
        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));
        
        return(id);
      }
    
    
    unsigned int
    insert_growth_params(sqlite3* db, transaction_manager& mgr, const growth_params& data,
                         const sqlite3_policy& policy)
      {
        assert(db != nullptr);
        
        // get number of rows in table; this will be the identifier for the new data set
        unsigned int new_id = count(db, tokenization_table<growth_params_token>(policy));
        
        std::ostringstream insert_stmt;
        insert_stmt
          << "INSERT INTO " << tokenization_table<growth_params_token>(policy) << " VALUES (@id, @abs, @rel);";
        
        // prepare SQL statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));
        
        // bind values to the parameters
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@id"), new_id));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@abs"), data.get_abserr()));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@rel"), data.get_relerr()));
        
        // perform insertion
        check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_GROWTH_PARAMS_FAIL, SQLITE_DONE);
        
        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));
        
        return(new_id);
      }
    
    
  }   // namespace sqlite3_operations
