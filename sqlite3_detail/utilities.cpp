//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include <iostream>
#include <sstream>
#include <assert.h>

#include "utilities.h"

#include "exceptions.h"
#include "localizations/en_GB/en_GB.h"

namespace sqlite3_operations
  {

    // error-check an exec statement
    void exec(sqlite3* db, const std::string& stmt, const std::string& err)
      {
        char* errmsg;

        int status = sqlite3_exec(db, stmt.c_str(), nullptr, nullptr, &errmsg);

        if(status != SQLITE_OK)
          {
            std::ostringstream msg;
            msg << err << errmsg << ") [status=" << status << "]";
            throw runtime_exception(exception_type::sqlite3_error, msg.str());
          }
      }


    // error-check an exec statement
    void exec(sqlite3* db, const std::string& stmt)
      {
        char* errmsg;

        int status = sqlite3_exec(db, stmt.c_str(), nullptr, nullptr, &errmsg);

        if(status != SQLITE_OK)
          {
            std::ostringstream msg;
            msg << ERROR_SQLITE3 << " " << errmsg << " [status=" << status << "]";
            throw runtime_exception(exception_type::sqlite3_error, msg.str());
          }
      }


    // error-check a non-exec statement
    void check_stmt(sqlite3* db, int status, const std::string& err, int check_code)
      {
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
        if(status != check_code)
          {
            std::ostringstream msg;
            msg << ERROR_SQLITE3 << " " << sqlite3_errmsg(db) << " [status=" << status << "]";
            std::cerr << msg.str() << '\n';
            assert(false);
            throw runtime_exception(exception_type::sqlite3_error, msg.str());
          }
      }

  }   // namespace sqlite3_operations
