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


// TOKENIZATION

// TODO: consider simplifying this code via traits and templates


std::unique_ptr<FRW_model_token> data_manager::tokenize(const FRW_model& obj)
  {
    // open a new transaction on the database
    std::shared_ptr<transaction_manager> transaction = this->open_transaction();
    
    // lookup id for this model, or generate one if it does not already exist
    std::unique_ptr<FRW_model_token> id = this->tokenize(*transaction, obj);
    
    // commit the transaction
    transaction->commit();
    
    return std::move(id);
  }


std::unique_ptr<FRW_model_token> data_manager::tokenize(transaction_manager& mgr, const FRW_model& obj)
  {
    // lookup id for this model, or generate one if it does not already exist
    unsigned int id = this->lookup_or_insert(mgr, obj);
    return std::make_unique<FRW_model_token>(id);
  }


std::unique_ptr<z_token> data_manager::tokenize(double z)
  {
    // open a new transaction on the database
    std::shared_ptr<transaction_manager> transaction = this->open_transaction();
    
    // lookup id for this redshift, or generate one if it does not already exist
    std::unique_ptr<z_token> id = this->tokenize(*transaction, z);
    
    // commit the transaction
    transaction->commit();
    
    return std::move(id);
  }


std::unique_ptr<z_token> data_manager::tokenize(transaction_manager& mgr, double z)
  {
    // lookup id for this redshift, or generate one if it does not already exist
    unsigned int id = this->lookup_or_insert(mgr, z);
    return std::make_unique<z_token>(id);
  }


std::unique_ptr<linear_Pk_token>
data_manager::tokenize(const FRW_model_token& model, const linear_Pk& Pk_lin)
  {
    // open a new transaction on the database
    std::shared_ptr<transaction_manager> transaction = this->open_transaction();
    
    // lookup id for this power spectrum, or generate one if it doesn't already exist
    std::unique_ptr<linear_Pk_token> id = this->tokenize(*transaction, model, Pk_lin);
    
    // commit the transaction
    transaction->commit();
    
    return std::move(id);
  }


std::unique_ptr<linear_Pk_token>
data_manager::tokenize(transaction_manager& mgr, const FRW_model_token& model, const linear_Pk& Pk_lin)
  {
    // lookup id for this power spectrum, or generate one if it doesn't already exist
    unsigned int id = this->lookup_or_insert(mgr, model, Pk_lin);
    return std::make_unique<linear_Pk_token>(id);
  }
