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


template <>
std::unique_ptr<oneloop_growth> data_manager::find<oneloop_growth>(transaction_manager& mgr, const FRW_model_token& model, const z_database& z_db)
  {
    return sqlite3_operations::find(this->handle, mgr, this->policy, model, z_db);
  }


template <>
std::unique_ptr<initial_filtered_Pk> data_manager::find<initial_filtered_Pk>(transaction_manager& mgr, const linear_Pk_token& token, const k_database& k_db)
  {
    return sqlite3_operations::find<initial_filtered_Pk>(this->handle, mgr, this->policy, token, k_db);
  }


template <>
std::unique_ptr<final_filtered_Pk> data_manager::find<final_filtered_Pk>(transaction_manager& mgr, const linear_Pk_token& token, const k_database& k_db)
  {
    return sqlite3_operations::find<final_filtered_Pk>(this->handle, mgr, this->policy, token, k_db);
  }


template <>
std::unique_ptr<loop_integral>
data_manager::find<loop_integral>(transaction_manager& mgr, const FRW_model_token& model, const k_token& k, const linear_Pk_token& Pk,
                                  const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff)
  {
    return sqlite3_operations::find(this->handle, mgr, this->policy, model, k, Pk, IR_cutoff, UV_cutoff);
  }


template <>
std::unique_ptr<oneloop_Pk>
data_manager::find<oneloop_Pk>(transaction_manager& mgr, const FRW_model_token& model, const k_token& k, const z_token& z,
                               const linear_Pk_token& Pk, const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff)
  {
    return sqlite3_operations::find(this->handle, mgr, this->policy, model, k, z, Pk, IR_cutoff, UV_cutoff);
  }


template <>
std::unique_ptr<Matsubara_XY>
data_manager::find(transaction_manager& mgr, const FRW_model_token& model, const linear_Pk_token& Pk,
                   const IR_resum_token& IR_resum)
  {
    return sqlite3_operations::find(this->handle, mgr, this->policy, model, Pk, IR_resum);
  }
