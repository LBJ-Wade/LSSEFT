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

constexpr auto ERROR_SQLITE3_MULTIPLE_PK_LINEAR        = "multiple linear power spectra with matching values";
constexpr auto ERROR_SQLITE3_PK_LINEAR                 = "linear power spectrum";
constexpr auto ERROR_SQLITE3_PK_LINEAR_WRONG_MODEL     = "was originally tagged for use with model id";
constexpr auto ERROR_SQLITE3_PK_LINEAR_WRONG_MD5       = "has changed MD5 hash";
constexpr auto ERROR_SQLITE3_INSERT_PK_LINEAR_FAIL     = "failed to insert linear Pk record [backend code=";

constexpr auto ERROR_SQLITE3_INSERT_TRANSFER_FAIL      = "failed to insert transfer function record [backend code=";
constexpr auto ERROR_SQLITE3_INSERT_GROWTH_G_FAIL      = "failed to insert one-loop growth g-factor record [backend code=";
constexpr auto ERROR_SQLITE3_INSERT_GROWTH_F_FAIL      = "failed to insert one-loop growth g-factor record [backend code=";
constexpr auto ERROR_SQLITE3_INSERT_LOOP_MOMENTUM_FAIL = "failed to insert one-loop momentum integral record [backend code=";
constexpr auto ERROR_SQLITE3_INSERT_ONELOOP_PK_FAIL    = "failed to insert one-loop P(k) record [backend code=";
constexpr auto ERROR_SQLITE3_INSERT_ONELOOP_RSD_PK_FAIL = "failed to insert one-loop RSD P(k) record [backend code=";
constexpr auto ERROR_SQLITE3_INSERT_MULTIPOLE_PK_FAIL  = "failed to insert multipole P(k) record [backend code=";
constexpr auto ERROR_SQLITE3_INSERT_MATSUBARA_A_FAIL   = "failed to insert Matsubara-A record [backend code=";

constexpr auto ERROR_SQLITE3_FG_GROWTH_TABLE_READ_FAIL = "failed to read from g- and f-factor growth tables [backend code=";
constexpr auto ERROR_SQLITE3_FG_GROWTH_MISREAD         = "read unexpected number of results from g- and f-factor growth table";
constexpr auto ERROR_SQLITE3_READ_LOOP_MOMENTUM_FAIL   = "failed to read from loop momentum table [backend code=";
constexpr auto ERROR_SQLITE3_LOOP_MOMENTUM_MISREAD     = "read unexpected number of results from loop momentum table [backend code=";
constexpr auto ERROR_SQLITE3_READ_PK_FAIL              = "failed to read from the delta-delta P(k) table [backend code=";
constexpr auto ERROR_SQLITE3_READ_PK_MISREAD           = "read unexpected number of results from delta-delta P(k) table [backend code=";
constexpr auto ERROR_SQLITE3_READ_RSD_PK_FAIL          = "failed to read from a delta-delta RSD P(k) table [backend code=";
constexpr auto ERROR_SQLITE3_READ_RSD_PK_MISREAD       = "read unexpected number of results from a delta-delta RSD P(k) table [backend code=";
constexpr auto ERROR_SQLITE3_READ_MATSUBARA_A_FAIL     = "failed to read from Matsubara-A table [backend code=";
constexpr auto ERROR_SQLITE3_MATSUBARA_A_MISREAD       = "read unexpected number of results from Matsubara-A table [backend code=";

#endif //LSSEFT_SQLITE3_EN_GB_H
