//
// Created by David Seery on 09/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
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
std::unique_ptr<wiggle_Pk> data_manager::find<wiggle_Pk>(transaction_manager& mgr, const linear_Pk_token& token, const k_database& k_db)
  {
    return sqlite3_operations::find(this->handle, mgr, this->policy, token, k_db);
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
