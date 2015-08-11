//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "transaction_manager.h"

#include "exceptions.h"
#include "localizations/en_GB/en_GB.h"


transaction_manager::transaction_manager(transaction_manager::open_handler& o, transaction_manager::commit_handler& c,
                                         transaction_manager::rollback_handler& r, transaction_manager::release_handler& l)
  : rolled_back(false),
    committed(false),
    do_open(o),
    do_commit(c),
    do_rollback(r),
    do_release(l)
  {
    // call open handler
    do_open();
  }


transaction_manager::~transaction_manager()
  {
    // rollback the transaction if it was not committed
    if(!this->committed) this->rollback();
  }


void transaction_manager::commit()
  {
    if(this->committed) return;
    if(this->rolled_back) throw runtime_exception(exception_type::transaction_error, ERROR_COMMIT_AFTER_ROLLBACK);

    this->do_commit();
    this->committed = true;

    this->do_release();
  }


void transaction_manager::rollback()
  {
    if(this->rolled_back) return;
    if(this->committed) throw runtime_exception(exception_type::transaction_error, ERROR_ROLLBACK_AFTER_COMMIT);

    this->do_rollback();
    this->rolled_back = true;

    this->do_release();
  }
