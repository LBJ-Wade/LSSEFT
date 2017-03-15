//
// Created by David Seery on 11/08/2015.
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


#include <iostream>
#include <sstream>
#include <assert.h>

#include "utilities.h"

#include "exceptions.h"
#include "localizations/messages.h"

#include "boost/optional.hpp"


namespace sqlite3_operations
  {

    // error-check an exec statement
    void exec(sqlite3* db, const std::string& stmt, const std::string& err)
      {
        assert(db != nullptr);

        char* errmsg;
        int status = sqlite3_exec(db, stmt.c_str(), nullptr, nullptr, &errmsg);

        if(status != SQLITE_OK)
          {
            std::ostringstream msg;
            msg << err << errmsg << ") [status=" << status << ", SQL=\"" << stmt << "\"]";
            throw runtime_exception(exception_type::sqlite3_error, msg.str());
          }
      }


    // error-check an exec statement
    void exec(sqlite3* db, const std::string& stmt)
      {
        assert(db != nullptr);

        char* errmsg;
        int status = sqlite3_exec(db, stmt.c_str(), nullptr, nullptr, &errmsg);

        if(status != SQLITE_OK)
          {
            std::ostringstream msg;
            msg << ERROR_SQLITE3 << " " << errmsg << " [status=" << status << ", SQL=\"" << stmt << "\"]";
            throw runtime_exception(exception_type::sqlite3_error, msg.str());
          }
      }


    // error-check a non-exec statement
    void check_stmt(sqlite3* db, int status, const std::string& err, int check_code)
      {
        assert(db != nullptr);

        if(status != check_code)
          {
            std::ostringstream msg;
            msg << err << sqlite3_errmsg(db) << ") [status=" << status << "]";
            throw runtime_exception(exception_type::sqlite3_error, msg.str());
          }
      }


    // error-check a non-exec statement
    void check_stmt(sqlite3* db, int status, int check_code)
      {
        assert(db != nullptr);

        if(status != check_code)
          {
            std::ostringstream msg;
            msg << ERROR_SQLITE3 << " " << sqlite3_errmsg(db) << " [status=" << status << "]";
            throw runtime_exception(exception_type::sqlite3_error, msg.str());
          }
      }


    unsigned int count(sqlite3* db, const std::string& table)
      {
        assert(db != nullptr);

        std::ostringstream count_stmt;
        count_stmt
          << "SELECT COUNT(*) FROM " << table << ";";

        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, count_stmt.str().c_str(), count_stmt.str().length()+1, &stmt, nullptr));

        boost::optional<unsigned int> rows;
        int status = 0;
        while((status = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(status == SQLITE_ROW)
              {
                if(rows) throw runtime_exception(exception_type::sqlite3_error, ERROR_SQLITE3_MULTIPLE_COUNT_ROWS);
                rows = static_cast<unsigned int>(sqlite3_column_int(stmt, 0));
              }
          }

        if(!rows) throw runtime_exception(exception_type::sqlite3_error, ERROR_SQLITE3_NO_COUNT_ROWS);

        return(*rows);
      }

  }   // namespace sqlite3_operations
