//
// Created by David Seery on 14/11/2017.
// --@@
// Copyright (c) 2017 University of Sussex. All rights reserved.
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


#include <sstream>

#include "pipeline_id.h"
#include "utilities.h"

#include "exceptions.h"
#include "localizations/messages.h"


namespace sqlite3_operations
  {

    void write_pipeline_id(sqlite3* db, const sqlite3_policy& policy, std::string tag)
      {
        std::ostringstream insert_stmt;
        insert_stmt << "INSERT INTO " << policy.pipeline_id_table() << " VALUES(@pipeline_id);";

        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));

        // bind parameter values
        check_stmt(db, sqlite3_bind_text(stmt, sqlite3_bind_parameter_index(stmt, "@pipeline_id"), tag.c_str(), tag.length(), SQLITE_STATIC));

        // perform insertion
        check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_PIPELINE_ID_FAIL, SQLITE_DONE);

        // clear bindings and release
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));
      }


    std::string read_pipeline_id(sqlite3* db, const sqlite3_policy& policy)
      {
        std::ostringstream read_stmt;
        read_stmt << "SELECT pipeline_id FROM " << policy.pipeline_id_table() << ";" << '\n';

        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, read_stmt.str().c_str(), read_stmt.str().length()+1, &stmt, nullptr));

        // perform read
        int result = 0;
        std::string pipeline;
        unsigned int count = 0;
        while((result = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(result == SQLITE_ROW)
              {
                std::string x(reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0)));
                pipeline = x;
                ++count;
              }
            else
              {
                check_stmt(db, sqlite3_clear_bindings(stmt));
                check_stmt(db, sqlite3_finalize(stmt));

                throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_READ_PIPELINE_ID_FAIL);
              }
          }

        // clear bindings and release
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));

        if(count != 1) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_PIPELINE_ID_MISREAD);

        return pipeline;
      }

  }
