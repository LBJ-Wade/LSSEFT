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

#ifndef LSSEFT_SQLITE3_FRW_MODEL_H
#define LSSEFT_SQLITE3_FRW_MODEL_H


#include <memory>

#include "cosmology/FRW_model.h"

#include "database/transaction_manager.h"
#include "sqlite3_policy.h"

#include "boost/optional.hpp"

#include "sqlite3.h"


namespace sqlite3_operations
  {

    //! lookup id for an FRW model
    boost::optional<unsigned int> lookup_FRW_model(sqlite3* db, transaction_manager& mgr,
                                                   const FRW_model& obj, const sqlite3_policy& policy, double tol);


    //! insert an FRW model
    unsigned int insert_FRW_model(sqlite3* db, transaction_manager& mgr, const FRW_model& obj,
                                  const sqlite3_policy& policy);

  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_FRW_MODEL_H
