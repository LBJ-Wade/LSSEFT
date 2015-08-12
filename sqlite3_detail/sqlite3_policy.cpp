//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "sqlite3_policy.h"

#include "sqlite3_defaults.h"


sqlite3_policy::sqlite3_policy()
  : FRW_model(SQLITE3_DEFAULT_FRW_TABLE_NAME),
    time_config(SQLITE3_DEFAULT_TIME_CONFIGURATION_TABLE),
    wavenumber_config(SQLITE3_DEFAULT_WAVENUMBER_CONFIGURATION_TABLE)
  {
  }
