//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_EN_GB_H
#define LSSEFT_SQLITE3_EN_GB_H


#define ERROR_SQLITE3                        "sqlite3 error:"

#define ERROR_SQLITE3_MULTIPLE_COUNT_ROWS    "multiple rows returned from SQL COUNT query"
#define ERROR_SQLITE3_NO_COUNT_ROWS          "no rows returned from SQL COUNT query"

#define ERROR_SQLITE3_MULTIPLE_FRW_MODELS    "multiple FRW models with matching parameters"
#define ERROR_SQLITE3_INSERT_FRW_MODEL_FAIL  "failed to insert FRW model record [backend code="

#define ERROR_SQLITE3_MULTIPLE_REDSHIFTS     "multiple redshifts with matching values"
#define ERROR_SQLITE3_INSERT_REDSHIFT_FAIL   "failed to insert redshift record [backend code="

#define ERROR_SQLITE3_MULTIPLE_WAVENUMBERS   "multiple wavenumbers with matching values"
#define ERROR_SQLITE3_INSERT_WAVENUMBER_FAIL "failed to insert wavenumber record [backend code="

#define ERROR_SQLITE3_TEMPORARY_REDSHIFT     "failed to insert redshift record in temporary table [backend code="
#define ERROR_SQLITE3_TEMPORARY_WAVENUMBER   "failed to insert wavenumber record in temporary table [backend code="

#define ERROR_SQLITE3_INSERT_TRANSFER_FAIL   "failed to insert transfer function record [backend code="
#define ERROR_SQLITE3_INSERT_ONELOOP_FAIL    "failed to insert growth-factor and one-loop kernel function record [backend code="


#endif //LSSEFT_SQLITE3_EN_GB_H
