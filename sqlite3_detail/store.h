//
// Created by David Seery on 17/08/2015.
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

#ifndef LSSEFT_SQLITE3_STORE_H
#define LSSEFT_SQLITE3_STORE_H


#include "database/transaction_manager.h"

#include "cosmology/concepts/transfer_function.h"
#include "cosmology/concepts/filtered_Pk_value.h"
#include "cosmology/concepts/oneloop_growth.h"
#include "cosmology/concepts/loop_integral.h"
#include "cosmology/concepts/oneloop_Pk.h"
#include "cosmology/concepts/multipole_Pk.h"
#include "cosmology/concepts/Matsubara_XY.h"

#include "sqlite3_policy.h"

#include "sqlite3.h"


namespace sqlite3_operations
  {

    //! store a transfer function sample
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const transfer_function& sample);
    
    //! store a filtered power spectrum sample
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token&, const filtered_Pk_value& sample);

    //! store a one-loop growth factor sample
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const oneloop_growth& sample);

    //! store a loop momentum integral sample
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const loop_integral& sample);
    
    //! store Matsubara X & Y coefficients
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const Matsubara_XY& sample);
    
    //! store a one-loop Pk sample
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const std::list<oneloop_Pk_set>& sample);

    //! store a multipole Pk sample
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const multipole_Pk_set& sample);

    //! store a counterterm sample
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const multipole_counterterm_set& sample);

  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_STORE_H
