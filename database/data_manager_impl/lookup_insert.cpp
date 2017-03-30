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


// LOOKUP AND INSERT


unsigned int data_manager::lookup_or_insert(transaction_manager& mgr, const FRW_model& obj)
  {
    boost::optional<unsigned int> id = sqlite3_operations::lookup_FRW_model(this->handle, mgr, obj, this->policy, this->FRW_model_tol);
    if(id) return *id;
    
    return sqlite3_operations::insert_FRW_model(this->handle, mgr, obj, this->policy);
  }


unsigned int data_manager::lookup_or_insert(transaction_manager& mgr, double z)
  {
    boost::optional<unsigned int> id = sqlite3_operations::lookup_redshift(this->handle, mgr, z, this->policy, this->z_tol);
    if(id) return *id;
    
    return sqlite3_operations::insert_redshift(this->handle, mgr, z, this->policy);
  }


unsigned int data_manager::lookup_or_insert(transaction_manager& mgr, const Pk_filter_params& data)
  {
    boost::optional<unsigned int> id = sqlite3_operations::lookup_filter_params(this->handle, mgr, data, this->policy, this->filter_tol);
    if(id) return *id;
    
    return sqlite3_operations::insert_filter_params(this->handle, mgr, data, this->policy);
  }


unsigned int data_manager::lookup_or_insert(transaction_manager& mgr, const loop_integral_params& data)
  {
    boost::optional<unsigned int> id = sqlite3_operations::lookup_oneloop_params(this->handle, mgr, data, this->policy, this->oneloop_tol);
    if(id) return *id;
    
    return sqlite3_operations::insert_oneloop_params(this->handle, mgr, data, this->policy);
  }


unsigned int data_manager::lookup_or_insert(transaction_manager& mgr, const MatsubaraXY_params& data)
  {
    boost::optional<unsigned int> id = sqlite3_operations::lookup_MatsubaraXY_params(this->handle, mgr, data, this->policy, this->MatsubaraXY_tol);
    if(id) return *id;
    
    return sqlite3_operations::insert_MatsubaraXY_params(this->handle, mgr, data, this->policy);
  }


unsigned int data_manager::lookup_or_insert(transaction_manager& mgr, const growth_params& data)
  {
    boost::optional<unsigned int> id = sqlite3_operations::lookup_growth_params(this->handle, mgr, data, this->policy, this->growth_tol);
    if(id) return *id;
    
    return sqlite3_operations::insert_growth_params(this->handle, mgr, data, this->policy);
  }
