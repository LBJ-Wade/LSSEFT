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


// TRANSACTIONS


std::shared_ptr<transaction_manager> data_manager::open_transaction()
  {
    // check whether a transaction is already in progress; if so, raise an exception
    std::shared_ptr<transaction_manager> check = this->current_transaction.lock();
    if(check) throw runtime_exception(exception_type::transaction_error, ERROR_TRANSACTION_IN_PROGRESS);
    check.reset();
    
    // create a new transaction manager
    transaction_manager::open_handler     do_open     = std::bind(&data_manager::begin_transaction, this);
    transaction_manager::commit_handler   do_commit   = std::bind(&data_manager::commit_transaction, this);
    transaction_manager::rollback_handler do_rollback = std::bind(&data_manager::rollback_transaction, this);
    transaction_manager::release_handler  do_release  = std::bind(&data_manager::release_transaction, this);
    
    std::shared_ptr<transaction_manager> transaction = std::make_shared<transaction_manager>(do_open, do_commit, do_rollback, do_release);
    
    // record this transaction
    this->current_transaction = transaction;
    
    return(transaction);
  }


void data_manager::begin_transaction()
  {
    assert(this->handle != nullptr);
    sqlite3_operations::exec(this->handle, "BEGIN TRANSACTION");
  }


void data_manager::commit_transaction()
  {
    assert(this->handle != nullptr);
    sqlite3_operations::exec(this->handle, "COMMIT");
  }


void data_manager::rollback_transaction()
  {
    assert(this->handle != nullptr);
    sqlite3_operations::exec(this->handle, "ROLLBACK");
  }


void data_manager::release_transaction()
  {
    assert(this->handle != nullptr);
    
    // check whether a transaction is already in progress
    std::shared_ptr<transaction_manager> check = this->current_transaction.lock();
    if(!check) throw runtime_exception(exception_type::transaction_error, ERROR_NO_TRANSACTION_IN_PROGRESS);
    check.reset();
    
    this->current_transaction.reset();
  }
