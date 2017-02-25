//
// Created by David Seery on 11/08/2015.
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

#include "transaction_manager.h"

#include "exceptions.h"
#include "localizations/messages.h"


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
