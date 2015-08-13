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
    delta_m(SQLITE3_DEFAULT_DELTA_M_TABLE),
    delta_r(SQLITE3_DEFAULT_DELTA_R_TABLE),
    theta_m(SQLITE3_DEFAULT_THETA_M_TABLE),
    theta_r(SQLITE3_DEFAULT_THETA_R_TABLE),
    Phi(SQLITE3_DEFAULT_PHI_TABLE),
    oneloop(SQLITE3_DEFAULT_ONELOOP_TABLE)
  {
  }
