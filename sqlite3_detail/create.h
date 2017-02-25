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

#ifndef LSSEFT_SQLITE3_CREATE_H
#define LSSEFT_SQLITE3_CREATE_H


#include <sstream>

#include "sqlite3_policy.h"

#include "sqlite3.h"


namespace sqlite3_operations
  {

    void create_tables(sqlite3* db, const sqlite3_policy& policy);

  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_CREATE_H
