//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_UTILITIES_H
#define LSSEFT_SQLITE3_UTILITIES_H


#include <string>

#include "sqlite3.h"


namespace sqlite3_operations
  {

    //! error-check an exec statement
    void exec(sqlite3* db, const std::string& stmt, const std::string& err);

    //! error-check an exec statement
    void exec(sqlite3* db, const std::string& stmt);

    // error-check a non-exec statement
    void check_stmt(sqlite3* db, int status, const std::string& err, int check_code=SQLITE_OK);

    // error-check a non-exec statement
    void check_stmt(sqlite3* db, int status, int check_code=SQLITE_OK);

  }




#endif //LSSEFT_SQLITE3_UTILITIES_H
