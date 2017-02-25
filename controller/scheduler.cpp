//
// Created by David Seery on 14/08/2015.
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


#include <assert.h>

#include "scheduler.h"


scheduler_detail::worker_data::worker_data()
  : initialized(false),
    active(true),
    assigned(false)
  {
  }


void scheduler_detail::worker_data::initialize(unsigned int n)
  {
    this->number = n;
    this->initialized = true;
  }


scheduler::scheduler(unsigned int N)
  : workers(N),
    waiting_for_initialization(N),
    active(0),
    unassigned(N)
  {
    worker_list.clear();
    worker_list.resize(workers);
  }


bool scheduler::is_assignable() const
  {
    return(this->unassigned > 0);
  }


void scheduler::initialize_worker(unsigned int n)
  {
    assert(this->waiting_for_initialization > 0);

    if(!this->worker_list[n].is_initialized())
      {
        this->worker_list[n].initialize(n);

        --this->waiting_for_initialization;
        ++this->active;
      }
    else
      {
        assert(false);
      }
  }


std::vector<unsigned int> scheduler::make_assignment()
  {
    std::vector<unsigned int> list;

    for(worker_list_type::const_iterator t = this->worker_list.begin(); t != this->worker_list.end(); ++t)
      {
        if(!t->is_assigned()) list.push_back(t->get_worker_number());
      }

    return(list);
  }


void scheduler::mark_assigned(unsigned int n)
  {
    assert(!this->worker_list[n].is_assigned());
    assert(this->unassigned > 0);

    this->worker_list[n].mark_assigned();
    --this->unassigned;
  }


void scheduler::mark_unassigned(unsigned int n)
  {
    assert(this->worker_list[n].is_assigned());

    this->worker_list[n].mark_unassigned();
    ++this->unassigned;
  }


void scheduler::mark_inactive(unsigned int n)
  {
    assert(this->worker_list[n].is_active());
    assert(!this->worker_list[n].is_assigned());
    assert(this->active > 0);

    this->worker_list[n].mark_inactive();
    this->active--;
  }
