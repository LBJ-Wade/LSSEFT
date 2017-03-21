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

#ifndef LSSEFT_SQLITE3_EN_GB_H
#define LSSEFT_SQLITE3_EN_GB_H


constexpr auto ERROR_SQLITE3                                         = "sqlite3 error:";

constexpr auto ERROR_SQLITE3_MULTIPLE_COUNT_ROWS                     = "multiple rows returned from SQL COUNT query";
constexpr auto ERROR_SQLITE3_NO_COUNT_ROWS                           = "no rows returned from SQL COUNT query";

constexpr auto ERROR_SQLITE3_MULTIPLE_FRW_MODELS                     = "multiple FRW models with matching parameters";
constexpr auto ERROR_SQLITE3_INSERT_FRW_MODEL_FAIL                   = "failed to insert FRW model record [backend code=";

constexpr auto ERROR_SQLITE3_MULTIPLE_REDSHIFTS                      = "multiple redshifts with matching values";
constexpr auto ERROR_SQLITE3_INSERT_REDSHIFT_FAIL                    = "failed to insert redshift record [backend code=";

constexpr auto ERROR_SQLITE3_MULTIPLE_WAVENUMBERS                    = "multiple wavenumbers with matching values";
constexpr auto ERROR_SQLITE3_INSERT_WAVENUMBER_FAIL                  = "failed to insert wavenumber record [backend code=";

constexpr auto ERROR_SQLITE3_MULTIPLE_FILTER_PARAMS                  = "multiple filter parameter sets with matching values";
constexpr auto ERROR_SQLITE3_INSERT_FILTER_PARAMS_FAIL               = "failed to insert filter parameter record [backend code=";

constexpr auto ERROR_SQLITE3_MULTIPLE_ONELOOP_PARAMS                 = "multiple oneloop parameter sets with matching values";
constexpr auto ERROR_SQLITE3_INSERT_ONELOOP_PARAMS_FAIL              = "failed to insert oneloop parameter record [backend code=";

constexpr auto ERROR_SQLITE3_MULTIPLE_MATSUBARAXY_PARAMS             = "multiple Matsubara X&Y parameter sets with matching values";
constexpr auto ERROR_SQLITE3_INSERT_MATSUBARAXY_PARAMS_FAIL          = "failed to insert Matsubara X&Y parameter record [backend code=";

constexpr auto ERROR_SQLITE3_MULTIPLE_GROWTH_PARAMS                  = "multiple growth-function parameter sets with matching values";
constexpr auto ERROR_SQLITE3_INSERT_GROWTH_PARAMS_FAIL               = "failed to insert growth-function parameter record [backend code=";

constexpr auto ERROR_SQLITE3_TEMPORARY_REDSHIFT                      = "failed to insert redshift record in temporary table [backend code=";
constexpr auto ERROR_SQLITE3_TEMPORARY_WAVENUMBER                    = "failed to insert wavenumber record in temporary table [backend code=";

constexpr auto ERROR_SQLITE3_MULTIPLE_PK_LINEAR                      = "multiple linear power spectra with matching values";
constexpr auto ERROR_SQLITE3_PK_LINEAR                               = "linear power spectrum";
constexpr auto ERROR_SQLITE3_PK_LINEAR_WRONG_MODEL                   = "was originally tagged for use with model id";
constexpr auto ERROR_SQLITE3_PK_LINEAR_WRONG_MD5                     = "has changed MD5 hash";
constexpr auto ERROR_SQLITE3_INSERT_PK_LINEAR_DATA_FAIL              = "failed to insert linear Pk data record [backend code=";
constexpr auto ERROR_SQLITE3_INSERT_PK_LINEAR_CONFIG_FAIL            = "failed to insert linear Pk configuration record [backend code=";

constexpr auto ERROR_SQLITE3_INSERT_TRANSFER_FAIL                    = "failed to insert transfer function record";
constexpr auto ERROR_SQLITE3_INSERT_GROWTH_D_FAIL                    = "failed to insert one-loop growth D-factor record";
constexpr auto ERROR_SQLITE3_INSERT_GROWTH_F_FAIL                    = "failed to insert one-loop growth f-factor record";
constexpr auto ERROR_SQLITE3_INSERT_LOOP_MOMENTUM_FAIL               = "failed to insert one-loop momentum integral record";
constexpr auto ERROR_SQLITE3_INSERT_ONELOOP_PK_FAIL                  = "failed to insert one-loop P(k) record";
constexpr auto ERROR_SQLITE3_INSERT_ONELOOP_RSD_PK_FAIL              = "failed to insert one-loop RSD P(k) record";
constexpr auto ERROR_SQLITE3_INSERT_RESUM_ONE_LOOP_PK_FAIL           = "failed to insert resummed one-loop P(k) record";
constexpr auto ERROR_SQLITE3_INSERT_MULTIPOLE_PK_FAIL                = "failed to insert multipole P(k) record";
constexpr auto ERROR_SQLITE3_INSERT_MATSUBARA_XY_FAIL                = "failed to insert Matsubara-XY record";

constexpr auto ERROR_SQLITE3_DF_GROWTH_TABLE_READ_FAIL               = "failed to read from D- and f-factor growth tables";
constexpr auto ERROR_SQLITE3_DF_GROWTH_MISREAD                       = "read unexpected number of results from D- and f-factor growth table";
constexpr auto ERROR_SQLITE3_READ_LOOP_MOMENTUM_FAIL                 = "failed to read from loop momentum table";
constexpr auto ERROR_SQLITE3_LOOP_MOMENTUM_MISREAD                   = "read unexpected number of results from loop momentum table";
constexpr auto ERROR_SQLITE3_READ_PK_FAIL                            = "failed to read from the delta-delta P(k) table";
constexpr auto ERROR_SQLITE3_READ_PK_MISREAD                         = "read unexpected number of results from delta-delta P(k) table";
constexpr auto ERROR_SQLITE3_READ_RSD_PK_FAIL                        = "failed to read from a delta-delta RSD P(k) table";
constexpr auto ERROR_SQLITE3_READ_RSD_PK_MISREAD                     = "read unexpected number of results from a delta-delta RSD P(k) table";
constexpr auto ERROR_SQLITE3_READ_MATSUBARA_XY_FAIL                  = "failed to read from Matsubara-XY table";
constexpr auto ERROR_SQLITE3_MATSUBARA_XY_MISREAD                    = "read unexpected number of results from Matsubara-XY table";
constexpr auto ERROR_SQLITE3_READ_FILTERED_PK_FAIL                   = "failed to read from filtered linear Pk table";
constexpr auto ERROR_SQLITE3_FILTERED_PK_MISREAD                     = "read unexpected number of results from filtered linear Pk table";

#endif //LSSEFT_SQLITE3_EN_GB_H
