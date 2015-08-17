//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "sqlite3_policy.h"

#include "sqlite3_defaults.h"


sqlite3_policy::sqlite3_policy()
  : FRW_model(SQLITE3_DEFAULT_FRW_TABLE_NAME),
    redshift_config(SQLITE3_DEFAULT_REDSHIFT_CONFIGURATION_TABLE),
    wavenumber_config(SQLITE3_DEFAULT_WAVENUMBER_CONFIGURATION_TABLE),
    transfer(SQLITE3_DEFAULT_TRANSFER_TABLE),
    oneloop(SQLITE3_DEFAULT_ONELOOP_TABLE),
    temp(SQLITE3_DEFAULT_TEMPORARY_TABLE)
  {
  }
