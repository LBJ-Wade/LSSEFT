//
// Created by David Seery on 12/08/2015.
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

#include "sqlite3_policy.h"

#include "sqlite3_defaults.h"


sqlite3_policy::sqlite3_policy()
  : FRW_model(SQLITE3_DEFAULT_FRW_TABLE_NAME),
    redshift_config(SQLITE3_DEFAULT_REDSHIFT_CONFIGURATION_TABLE),
    wavenumber_config(SQLITE3_DEFAULT_WAVENUMBER_CONFIGURATION_TABLE),
    IR_config(SQLITE3_DEFAULT_IR_CUTOFF_CONFIGURATION_TABLE),
    UV_config(SQLITE3_DEFAULT_UV_CUTOFF_CONFIGURATION_TABLE),
    IR_resum_config(SQLITE3_DEFAULT_IR_RESUM_CONFIGURATION_TABLE),
    Pk_linear_config(SQLITE3_DEFAULT_LINEAR_PK_CONFIGURATION_TABLE),
    Pk_linear(SQLITE3_DEFAULT_LINEAR_PK_DATA_TABLE),
    filter_config(SQLITE3_DEFAULT_FILTER_CONFIGURATION_TABLE),
    loop_integral_config(SQLITE3_DEFAULT_LOOP_INTEGRAL_CONFIGURATION_TABLE),
    MatsubaraXY_config(SQLITE3_DEFAULT_MATSUBARAXY_CONFIGURATION_TABLE),
    growth_config(SQLITE3_DEFAULT_GROWTH_CONFIGURATION_TABLE),
    transfer(SQLITE3_DEFAULT_TRANSFER_TABLE),
    growth_D_factor(SQLITE3_DEFAULT_GROWTH_D_TABLE),
    growth_f_factor(SQLITE3_DEFAULT_GROWTH_F_TABLE),
    Matsubara_XY(SQLITE3_DEFAULT_MATSUBARA_XY_TABLE),
    counterterms_c0(SQLITE3_DEFAULT_COUNTERTERMS_C0_TABLE),
    counterterms_c2(SQLITE3_DEFAULT_COUNTERTERMS_C2_TABLE),
    counterterms_c4(SQLITE3_DEFAULT_COUNTERTERMS_C4_TABLE),
    counterterms_c6(SQLITE3_DEFAULT_COUNTERTERMS_C6_TABLE),
    temp(SQLITE3_DEFAULT_TEMPORARY_TABLE)
  {
  }
