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
    
    
    void write_performance_pragmas(sqlite3* db, bool network_filesystem)
      {
        assert(db != nullptr);
        
        // SQLite performance choices:
        // http://blog.devart.com/increasing-sqlite-performance.html
        // https://wiki.mozilla.org/Performance/Avoid_SQLite_In_Your_Next_Firefox_Feature#Important_Pragmas
        // http://stackoverflow.com/questions/784173/what-are-the-performance-characteristics-of-sqlite-with-very-large-database-file
    
        // attempt to speed up insertions by disabling foreign key constraints
        // during initial writes
        char* errmsg;
        sqlite3_exec(db, "PRAGMA foreign_keys = OFF;", nullptr, nullptr, &errmsg);

        // if write-ahead log mode is disabled (as it must be if a network filing system is in play)
        // then put journal into truncate mode
        // otherwise, enable to write-ahead log
        if(network_filesystem)
          {
            sqlite3_exec(db, "PRAGMA journal_mode = TRUNCATE;", nullptr, nullptr, &errmsg);
          }
        else
          {
            sqlite3_exec(db, "PRAGMA journal_mode = WAL;", nullptr, nullptr, &errmsg);
          }
    
        // force temporary objects to be stored in memory, for speed
        sqlite3_exec(db, "PRAGMA temp_store = MEMORY;", nullptr, nullptr, &errmsg);
    
        // set SYNCHRONOUS mode to 'Normal' rather than 'Full'
        sqlite3_exec(db, "PRAGMA synchronous = NORMAL;", nullptr, nullptr, &errmsg);
    
        // CACHE_SIZE unlikely to make much difference except in windows
        sqlite3_exec(db, "PRAGMA cache_size = 10000;", nullptr, nullptr, &errmsg);
    
        // PAGE_SIZE unlikely to make much difference except in windows
        sqlite3_exec(db, "PRAGMA page_size = 4096;", nullptr, nullptr, &errmsg);
      }
    
    
    void default_pragmas(sqlite3* db)
      {
        assert(db != nullptr);
        
        // SQLite performance choices:
        // http://blog.devart.com/increasing-sqlite-performance.html
        // https://wiki.mozilla.org/Performance/Avoid_SQLite_In_Your_Next_Firefox_Feature#Important_Pragmas
        // http://stackoverflow.com/questions/784173/what-are-the-performance-characteristics-of-sqlite-with-very-large-database-file
    
        char* errmsg;
        sqlite3_exec(db, "PRAGMA foreign_keys = ON;", nullptr, nullptr, &errmsg);
    
        // don't change the journal mode
    
        // force temporary objects to be stored in memory, for speed
        sqlite3_exec(db, "PRAGMA temp_store = MEMORY;", nullptr, nullptr, &errmsg);
    
        // set SYNCHRONOUS mode to 'Normal' rather than 'Full'
        sqlite3_exec(db, "PRAGMA synchronous = NORMAL;", nullptr, nullptr, &errmsg);
    
        // CACHE_SIZE unlikely to make much difference except in windows
        sqlite3_exec(db, "PRAGMA cache_size = 10000;", nullptr, nullptr, &errmsg);
    
        // PAGE_SIZE unlikely to make much difference except in windows
        sqlite3_exec(db, "PRAGMA page_size = 4096;", nullptr, nullptr, &errmsg);
      }
    
    
    void create_index(sqlite3* db, const std::string& table, const std::string& column)
      {
        assert(db != nullptr);
        
        std::ostringstream index_stmt;
        index_stmt
          << "CREATE INDEX " << table << "_" << column << "_idx ON " << table << "(" << column << ");";
        exec(db, index_stmt.str());
      }
    
    
    void drop_index(sqlite3* db, const std::string& table, const std::string& column)
      {
        assert(db != nullptr);
    
        std::ostringstream index_stmt;
        index_stmt
          << "DROP INDEX IF EXISTS " << table << "_" << column << "_idx;";
        exec(db, index_stmt.str());
      }
    
    
    void analyze(sqlite3* db)
      {
        assert(db != nullptr);
        exec(db, "ANALYZE;");
        exec(db, "ANALYZE sqlite_master;");
      }
    
    
    void create_index(sqlite3* db, const std::string& table, std::initializer_list<std::string> list)
      {
        for(const std::string& col : list)
          {
            create_index(db, table, col);
          }
      }
    
    
    void drop_index(sqlite3* db, const std::string& table, std::initializer_list<std::string> list)
      {
        for(const std::string& col : list)
          {
            drop_index(db, table, col);
          }
      }
    
    
  }   // namespace sqlite3_operations
