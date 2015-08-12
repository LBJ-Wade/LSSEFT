//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_CREATE_H
#define LSSEFT_SQLITE3_CREATE_H


#include <sstream>

#include "sqlite3_policy.h"

#include "sqlite3.h"


namespace sqlite3_operations
  {

    void create_tables(sqlite3* db, const sqlite3_policy& policy);

  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_CREATE_H
