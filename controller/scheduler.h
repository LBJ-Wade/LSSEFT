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

#ifndef LSSEFT_SCHEDULER_H
#define LSSEFT_SCHEDULER_H


#include <vector>


namespace scheduler_detail
  {

    class worker_data
      {

        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! constructor
        worker_data();

        //! destructor is defaults
        ~worker_data() = default;


        // PROPERTIES

      public:

        //! is this worker initialized?
        bool is_initialized() const { return(this->initialized); }

        //! is this worker active?
        bool is_active() const { return(this->active); }

        //! is this worker assigned
        bool is_assigned() const { return(this->assigned); }

        //! get worker number
        unsigned int get_worker_number() const { if(this->initialized) return(this->number); else return(0); }


        // MANAGEMENT

      public:

        //! set initialization
        void initialize(unsigned int n);

        //! set active status
        void mark_inactive() { this->active = false; }

        //! set assigned status
        void mark_assigned() { this->assigned = true; }

        //! set unassigned status
        void mark_unassigned() { this->assigned = false; }


        // INTERNAL DATA

      private:

        //! worker number
        unsigned int number;

        //! has this worker been initialized?
        bool initialized;

        //! is this worker active?
        bool active;

        //! is this worker assigned
        bool assigned;

      };

  }

class scheduler
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! construct new scheduler for N workers
    scheduler(unsigned int N);

    //! destructor is default
    ~scheduler() = default;


    // QUERY STATUS

  public:

    //! ready for scheduling to proceed?
    bool is_ready() const { return(this->waiting_for_initialization == 0); }

    //! can work be assigned?
    bool is_assignable() const;

    //! are all workers inactive?
    bool all_inactive() const { return(this->active == 0); }


    // WORKER MANAGEMENT

  public:

    //! initialize a worker
    void initialize_worker(unsigned int n);

    //! mark a worker assigned
    void mark_assigned(unsigned int n);

    //! mark a worker unassigned
    void mark_unassigned(unsigned int n);

    //! mark a worker inactive
    void mark_inactive(unsigned int n);


    // WORK ASSIGNMENTS

  public:

    //! construct list workers requiring assignment
    std::vector<unsigned int> make_assignment();


    // INTERNAL DATA

  private:


    // QUICK-ACCESS DATA

    //! number of worker processes
    unsigned int workers;

    //! number of workers waiting for initialization
    unsigned int waiting_for_initialization;

    //! number of active workers
    unsigned int active;

    //! number of unassigned workers
    unsigned int unassigned;

    // WORKER MANAGEMENT

    //! type alias for list of workers
    typedef std::vector< scheduler_detail::worker_data > worker_list_type;

    //! list of workers
    worker_list_type worker_list;

  };


#endif //LSSEFT_SCHEDULER_H
