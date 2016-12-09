//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include <iostream>
#include <sstream>
#include <assert.h>

#include <set>
#include <unordered_set>

#include "database/data_manager.h"
#include "database/data_manager_impl/types.h"

#include "sqlite3_detail/utilities.h"
#include "sqlite3_detail/operations.h"

#include "utilities/formatter.h"

#include "defaults.h"

#include "boost/timer/timer.hpp"


data_manager::data_manager(const boost::filesystem::path& c)
  : container(c),
    handle(nullptr),   // try to catch handle-not-initialized errors
    policy(),
    FRW_model_tol(LSSEFT_DEFAULT_FRW_MODEL_PARAMETER_TOLERANCE),
    z_tol(LSSEFT_DEFAULT_REDSHIFT_CONFIGURATION_TOLERANCE),
    k_tol(LSSEFT_DEFAULT_WAVENUMBER_CONFIGURATION_TOLERANCE)
  {
    // check whether container already exists
    if(boost::filesystem::exists(container))
      {
        if(boost::filesystem::is_regular_file(container))
          {
            if(sqlite3_open_v2(container.string().c_str(), &handle, SQLITE_OPEN_READWRITE, nullptr) != SQLITE_OK)
              {
                std::ostringstream msg;
                msg << ERROR_DATABASE_SQLITE_OPEN_FAILED << " " << container;
                throw runtime_exception(exception_type::database_error, msg.str());
              }

            return;
          }
        else
          {
            std::ostringstream msg;
            msg << ERROR_DATABASE_IS_NOT_FILE_A << " " << container << " " << ERROR_DATABASE_IS_NOT_FILE_B;
            throw runtime_exception(exception_type::database_error, msg.str());
          }
      }

    // if we get to here, the container does not already exist so we should create it
    if(sqlite3_open_v2(container.string().c_str(), &handle, SQLITE_OPEN_CREATE | SQLITE_OPEN_READWRITE, nullptr) != SQLITE_OK)
      {
        std::ostringstream msg;
        msg << ERROR_DATABASE_SQLITE_CREATE_FAILED << " " << container;
        throw runtime_exception(exception_type::database_error, msg.str());
      }

    // set up tables
    sqlite3_operations::create_tables(handle, policy);
  }


data_manager::~data_manager()
  {
    assert(this->handle != nullptr);

    // perform routine maintenance on the container
    sqlite3_operations::exec(this->handle, "VACUUM;");

    sqlite3_close(this->handle);
  }
