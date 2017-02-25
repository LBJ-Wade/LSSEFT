//
// Created by David Seery on 13/08/2015.
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

#ifndef LSSEFT_SQLITE3_TEMPORARY_TABLES_H
#define LSSEFT_SQLITE3_TEMPORARY_TABLES_H


#include <memory>
#include <string>

#include "sqlite3_policy.h"

#include "database/transaction_manager.h"
#include "database/tokens.h"
#include "database/z_database.h"
#include "database/k_database.h"

#include "sqlite3.h"


namespace sqlite3_operations
  {

    //! create temporary table of redshifts
    std::string z_table(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                        const z_database& z_db);

    //! create temporary table of wavenumbres
    std::string k_table(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                        const k_database& k_db);

    //! drop a temporary table
    void drop_temp(sqlite3* db, transaction_manager& mgr, const std::string& table);

  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_TEMPORARY_TABLES_H
