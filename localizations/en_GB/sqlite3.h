//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_EN_GB_H
#define LSSEFT_SQLITE3_EN_GB_H


constexpr auto ERROR_SQLITE3                           = "sqlite3 error:";

constexpr auto ERROR_SQLITE3_MULTIPLE_COUNT_ROWS       = "multiple rows returned from SQL COUNT query";
constexpr auto ERROR_SQLITE3_NO_COUNT_ROWS             = "no rows returned from SQL COUNT query";

constexpr auto ERROR_SQLITE3_MULTIPLE_FRW_MODELS       = "multiple FRW models with matching parameters";
constexpr auto ERROR_SQLITE3_INSERT_FRW_MODEL_FAIL     = "failed to insert FRW model record [backend code=";

constexpr auto ERROR_SQLITE3_MULTIPLE_REDSHIFTS        = "multiple redshifts with matching values";
constexpr auto ERROR_SQLITE3_INSERT_REDSHIFT_FAIL      = "failed to insert redshift record [backend code=";

constexpr auto ERROR_SQLITE3_MULTIPLE_WAVENUMBERS      = "multiple wavenumbers with matching values";
constexpr auto ERROR_SQLITE3_INSERT_WAVENUMBER_FAIL    = "failed to insert wavenumber record [backend code=";

constexpr auto ERROR_SQLITE3_TEMPORARY_REDSHIFT        = "failed to insert redshift record in temporary table [backend code=";
constexpr auto ERROR_SQLITE3_TEMPORARY_WAVENUMBER      = "failed to insert wavenumber record in temporary table [backend code=";

constexpr auto ERROR_SQLITE3_INSERT_TRANSFER_FAIL      = "failed to insert transfer function record [backend code=";
constexpr auto ERROR_SQLITE3_INSERT_GROWTH_G_FAIL      = "failed to insert one-loop growth g-factor record [backend code=";
constexpr auto ERROR_SQLITE3_INSERT_GROWTH_F_FAIL      = "failed to insert one-loop growth g-factor record [backend code=";
constexpr auto ERROR_SQLITE3_INSERT_GROWTH_META_FAIL   = "failed to insert one-loop growth metadata record [backend code=";
constexpr auto ERROR_SQLITE3_INSERT_LOOP_MOMENTUM_FAIL = "failed to insert one-loop momentum integral record [backend code=";

constexpr auto ERROR_SQLITE3_FG_GROWTH_TABLE_READ_FAIL = "failed to read from g- and f-factor growth tables [backend code=";
constexpr auto ERROR_SQLITE3_FG_GROWTH_META_READ_FAIL  = "failed to read g- and f-factor metadata [backend code=";
constexpr auto ERROR_SQLITE3_FG_GROWTH_MISREAD         = "read unexpected number of results from g- and f-factor growth table";
constexpr auto ERROR_SQLITE3_READ_LOOP_MOMENTUM_FAIL   = "failed to read from loop momentum table [backend code=";
constexpr auto ERROR_SQLITE3_LOOP_MOMENTUM_MISREAD     = "read unexpected number of results from loop momentum table [backend code=";

#endif //LSSEFT_SQLITE3_EN_GB_H
