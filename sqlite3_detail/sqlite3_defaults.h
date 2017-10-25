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

#ifndef LSSEFT_SQLITE3_DEFAULTS_H
#define LSSEFT_SQLITE3_DEFAULTS_H


constexpr auto SQLITE3_DEFAULT_FRW_TABLE_NAME                    = "models";
constexpr auto SQLITE3_DEFAULT_REDSHIFT_CONFIGURATION_TABLE      = "z_config";
constexpr auto SQLITE3_DEFAULT_WAVENUMBER_CONFIGURATION_TABLE    = "k_config";
constexpr auto SQLITE3_DEFAULT_IR_CUTOFF_CONFIGURATION_TABLE     = "IR_cutoff_config";
constexpr auto SQLITE3_DEFAULT_UV_CUTOFF_CONFIGURATION_TABLE     = "UV_cutoff_config";
constexpr auto SQLITE3_DEFAULT_IR_RESUM_CONFIGURATION_TABLE      = "IR_resum_config";
constexpr auto SQLITE3_DEFAULT_LINEAR_PK_CONFIGURATION_TABLE     = "Pk_linear_config";
constexpr auto SQLITE3_DEFAULT_LINEAR_PK_DATA_TABLE              = "Pk_linear";
constexpr auto SQLITE3_DEFAULT_FILTER_CONFIGURATION_TABLE        = "filter_config";
constexpr auto SQLITE3_DEFAULT_LOOP_INTEGRAL_CONFIGURATION_TABLE = "loop_config";
constexpr auto SQLITE3_DEFAULT_MATSUBARAXY_CONFIGURATION_TABLE   = "Matsubara_XY_config";
constexpr auto SQLITE3_DEFAULT_GROWTH_CONFIGURATION_TABLE        = "Df_factors_config";
constexpr auto SQLITE3_DEFAULT_TRANSFER_TABLE                    = "transfer";
constexpr auto SQLITE3_DEFAULT_GROWTH_D_TABLE                    = "D_factors";
constexpr auto SQLITE3_DEFAULT_GROWTH_F_TABLE                    = "f_factors";
constexpr auto SQLITE3_DEFAULT_MATSUBARA_XY_TABLE                = "Matsubara_XY";
constexpr auto SQLITE3_DEFAULT_COUNTERTERMS_C0_TABLE             = "counterterms_c0";
constexpr auto SQLITE3_DEFAULT_COUNTERTERMS_C2_TABLE             = "counterterms_c2";
constexpr auto SQLITE3_DEFAULT_COUNTERTERMS_C4_TABLE             = "counterterms_c4";
constexpr auto SQLITE3_DEFAULT_COUNTERTERMS_C6_TABLE             = "counterterms_c6";

constexpr auto SQLITE3_DEFAULT_TEMPORARY_TABLE                   = "temp";


#endif //LSSEFT_SQLITE3_DEFAULTS_H
