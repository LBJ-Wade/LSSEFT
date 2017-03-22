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

#ifndef LSSEFT_SQLITE3_UTILITIES_H
#define LSSEFT_SQLITE3_UTILITIES_H


#include <string>

#include "sqlite3.h"


namespace sqlite3_operations
  {

    // ERROR CHECKING

    //! error-check an exec statement
    void exec(sqlite3* db, const std::string& stmt, const std::string& err);

    //! error-check an exec statement
    void exec(sqlite3* db, const std::string& stmt);

    // error-check a non-exec statement
    void check_stmt(sqlite3* db, int status, const std::string& err, int check_code=SQLITE_OK);

    // error-check a non-exec statement
    void check_stmt(sqlite3* db, int status, int check_code=SQLITE_OK);


    // COUNTING FUNCTIONS

    //! count number of rows in a specified table
    unsigned int count(sqlite3* db, const std::string& table);
    
    
    // PRAGMA MANAGEMENT
    
    //! optimize SQLite performance settings for write performance
    void write_performance_pragmas(sqlite3* db, bool network_filesystem);
    
    //! relax SQLite settings
    void default_pragmas(sqlite3* db);
    
    
    // INDEX MANAGEMENT
    
    //! create a SQLite index
    void create_index(sqlite3* db, const std::string& table, const std::string& column);
    
    //! drop an SQLite index
    void drop_index(sqlite3* db, const std::string& table, const std::string& column);

    //! create a set of SQLite indices
    void create_index(sqlite3* db, const std::string& table, std::initializer_list<std::string> list);
    
    //! drop a set of SQLite indices
    void drop_index(sqlite3* db, const std::string& table, std::initializer_list<std::string> list);
    
    //! update SQLite's internal statistics
    void analyze(sqlite3* db);

  }




#endif //LSSEFT_SQLITE3_UTILITIES_H
