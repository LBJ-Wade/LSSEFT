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
