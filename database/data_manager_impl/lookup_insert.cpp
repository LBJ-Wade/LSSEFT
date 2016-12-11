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


// LOOKUP AND INSERT


unsigned int data_manager::lookup_or_insert(transaction_manager& mgr, const FRW_model& obj)
  {
    boost::optional<unsigned int> id = sqlite3_operations::lookup_FRW_model(this->handle, mgr, obj, this->policy, this->FRW_model_tol);
    if(id) return(*id);
    
    return sqlite3_operations::insert_FRW_model(this->handle, mgr, obj, this->policy);
  }


unsigned int data_manager::lookup_or_insert(transaction_manager& mgr, double z)
  {
    boost::optional<unsigned int> id = sqlite3_operations::lookup_redshift(this->handle, mgr, z, this->policy, this->z_tol);
    if(id) return(*id);
    
    return sqlite3_operations::insert_redshift(this->handle, mgr, z, this->policy);
  }


unsigned int data_manager::lookup_or_insert(transaction_manager& mgr, const FRW_model_token& model, const linear_Pk& Pk_lin)
  {
    boost::optional<unsigned int> id = sqlite3_operations::lookup_Pk_linear(this->handle, mgr, model, Pk_lin, this->policy);
    if(id) return(*id);
    
    return sqlite3_operations::insert_Pk_linear(this->handle, mgr, model, Pk_lin, this->policy);
  }
