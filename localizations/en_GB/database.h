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

#ifndef LSSEFT_DATABASE_EN_GB_H
#define LSSEFT_DATABASE_EN_GB_H


#define ERROR_DATABASE_IS_NOT_FILE_A        "cannot open specified database"
#define ERROR_DATABASE_IS_NOT_FILE_B        "because it is not a regular file"

#define ERROR_DATABASE_SQLITE_OPEN_FAILED   "failed to open specified database"
#define ERROR_DATABASE_SQLITE_CREATE_FAILED "failed to create specified database"

#define ERROR_TRANSACTION_IN_PROGRESS       "transaction already open"
#define ERROR_NO_TRANSACTION_IN_PROGRESS    "attempt to release transaction when none is currently open"

#define DATABASE_ATTACHED                   "attached to database"

#endif //LSSEFT_DATABASE_EN_GB_H
