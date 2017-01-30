//
// Created by David Seery on 09/12/2016.
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

#include <iostream>
#include <sstream>
#include <assert.h>

#include <set>
#include <unordered_set>

#include "database/data_manager.h"
#include "database/data_manager_impl/types.h"

#include "sqlite3_detail/utilities.h"
#include "sqlite3_detail/operations.h"

#include "utilities/formatter.h"

#include "defaults.h"

#include "boost/timer/timer.hpp"


std::unique_ptr<wiggle_Pk> data_manager::build_wiggle_Pk(const linear_Pk_token& token, const linear_Pk& Pk_lin)
  {
    // open a transaction on the database
    std::shared_ptr<transaction_manager> mgr = this->open_transaction();
    
    // extract database of wavenumber configurations from linear power spectrum container
    std::unique_ptr<k_database> k_db = this->build_k_db(*mgr, Pk_lin, FILTER_PK_DEFAULT_BOTTOM_CLEARANCE, FILTER_PK_DEFAULT_TOP_CLEARANCE);
    
    // extract wiggle_Pk container
    std::unique_ptr<wiggle_Pk> payload = this->find<wiggle_Pk>(*mgr, token, *k_db);
    
    // close transaction
    mgr->commit();
    
    return std::move(payload);
  }
