//
// Created by David Seery on 10/08/2015.
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

#ifndef LSSEFT_MASTER_CONTROLLER_H
#define LSSEFT_MASTER_CONTROLLER_H


#include <memory>

#include "argument_cache.h"
#include "local_environment.h"
#include "scheduler.h"

#include "database/data_manager.h"
#include "database/tokens.h"

#include "cosmology/types.h"
#include "cosmology/FRW_model.h"
#include "cosmology/oneloop_growth_integrator.h"
#include "cosmology/Pk_filter.h"
#include "cosmology/Matsubara_XY_calculator.h"
#include "cosmology/oneloop_growth_integrator.h"
#include "cosmology/concepts/range.h"
#include "cosmology/concepts/power_spectrum.h"

#include "error/error_handler.h"

#include "MPI_detail/mpi_traits.h"
#include "MPI_detail/mpi_payloads.h"

#include "boost/mpi.hpp"


class master_controller
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! construct a master controller delegate
    master_controller(boost::mpi::environment& me, boost::mpi::communicator& mw, argument_cache& ac);

    //! destructor is default
    ~master_controller() = default;


    // INTERFACE

  public:

    //! process command-line arguments
    void process_arguments(int argc, char* argv[]);

    //! execute
    void execute();


    // RANK TO WORKER NUMBER CONVERSIONS (worker number runs from 1 .. n-1, rank runs from 1 .. n, based on master process on rank 0)

  protected:

    //! Get worker number
    unsigned int worker_number() { return(static_cast<unsigned int>(this->mpi_world.rank()-1)); }

    //! Return MPI rank of this process
    unsigned int get_rank(void) const { return(static_cast<unsigned int>(this->mpi_world.rank())); }

    //! Map worker number to communicator rank
    unsigned int worker_rank(unsigned int worker_number) const { return(worker_number+1); }

    //! Map communicator rank to worker number
    unsigned int worker_number(unsigned int worker_rank) const { return(worker_rank-1); }


    // WORKER HANDLING

  protected:

    //! execute a job specified by a work list
    template <typename WorkItemList>
    void scatter(const FRW_model& model, const FRW_model_token& token, WorkItemList& work, data_manager& dmgr);

    //! store a payload returned by a worker
    template <typename WorkItem>
    void store_payload(const FRW_model_token& token, unsigned int source, data_manager& dmgr);

    //! terminate worker processes
    void terminate_workers();

    //! instruct workers to await new tasks
    std::unique_ptr<scheduler> set_up_workers(unsigned int tag);

    //! instruct workers that the current batch of tasks is finished
    void close_down_workers();


    // COMPUTE ONE-LOOP KERNELS

  protected:

    //! compute kernels at given redshifts
    void integrate_loop_growth(const FRW_model& model, const FRW_model_token& token, z_database& z_db, data_manager& dmgr,
                               const growth_params_token& params_tok, const growth_params& params);


    // INTERNAL DATA

  private:

    // MPI environment

    //! MPI environment object, inherited from parent task manager
    boost::mpi::environment& mpi_env;

    //! MPI communicator, inherited from parent task manager
    boost::mpi::communicator& mpi_world;


    // Caches

    //! argument cache, inherited from parent task manager
    argument_cache& arg_cache;

    //! local environment properties (constructed locally)
    local_environment local_env;


    // Functional blocks

    //! error handler
    error_handler err_handler;

  };


template <typename WorkItemList>
void master_controller::scatter(const FRW_model& model, const FRW_model_token& token, WorkItemList& work, data_manager& dmgr)
  {
    using WorkItem = typename WorkItemList::value_type;
    
    boost::timer::cpu_timer timer;              // total CPU time
    boost::timer::cpu_timer write_timer;        // time spent writing to the database
    write_timer.stop();
    
    if(this->mpi_world.size() == 1) throw runtime_exception(exception_type::runtime_error, ERROR_TOO_FEW_WORKERS);
    
    // ask data manager to prepare for new writes
    boost::timer::cpu_timer pre_timer;          // time spent doing preparation
    dmgr.setup_write(work);
    pre_timer.stop();
    
    // instruct slave processes to await transfer function tasks
    std::unique_ptr<scheduler> sch = this->set_up_workers(MPI_detail::work_item_traits<WorkItem>::new_task_message());
    
    bool sent_closedown = false;
    auto next_work_item = work.cbegin();
    
    while(!sch->all_inactive())
      {
        // check whether all work is exhausted
        if(next_work_item == work.cend() && !sent_closedown)
          {
            sent_closedown = true;
            this->close_down_workers();
          }
        
        // check whether any workers are waiting for assignments
        if(next_work_item != work.cend() && sch->is_assignable())
          {
            std::vector<unsigned int> unassigned_list = sch->make_assignment();
            std::vector<boost::mpi::request> requests;
            
            for(std::vector<unsigned int>::const_iterator t = unassigned_list.begin();
                next_work_item != work.cend() && t != unassigned_list.end(); ++t)
              {
                // assign next work item to this worker
                requests.push_back(this->mpi_world.isend(this->worker_rank(*t),
                                                         MPI_detail::work_item_traits<WorkItem>::new_item_message(),
                                                         MPI_detail::build_payload(model, next_work_item)));
                
                sch->mark_assigned(*t);
                ++next_work_item;
              }
            
            // wait for all messages to be received
            boost::mpi::wait_all(requests.begin(), requests.end());
          }
        
        // check whether any messages are waiting in the queue
        boost::optional<boost::mpi::status> stat = this->mpi_world.iprobe();
        
        while(stat) // consume messages until no more are available
          {
            switch(stat->tag())
              {
                case MPI_detail::MESSAGE_WORK_PRODUCT_READY:
                  {
                    write_timer.resume();
                    this->store_payload<WorkItem>(token, stat->source(), dmgr);
                    write_timer.stop();
                    sch->mark_unassigned(this->worker_number(stat->source()));
                    break;
                  }
                
                case MPI_detail::MESSAGE_END_OF_WORK_ACK:
                  {
                    this->mpi_world.recv(stat->source(), MPI_detail::MESSAGE_END_OF_WORK_ACK);
                    sch->mark_inactive(this->worker_number(stat->source()));
                    break;
                  }
                
                default:
                  {
                    assert(false);
                  }
              }
            
            stat = this->mpi_world.iprobe();
          }
      }
    
    boost::timer::cpu_timer post_timer;     // time spent tidying up the database after a write
    dmgr.finalize_write(work);
    post_timer.stop();
    
    timer.stop();
    std::ostringstream msg;
    msg << "completed work in time " << format_time(timer.elapsed().wall)
        << " ["
        << "database performance: prepare " << format_time(pre_timer.elapsed().wall) << ", "
        << "writes " << format_time(write_timer.elapsed().wall) << ", "
        << "cleanup " << format_time(post_timer.elapsed().wall)
        << "]";
    this->err_handler.info(msg.str());
  }


template <typename WorkItem>
void master_controller::store_payload(const FRW_model_token& token, unsigned int source, data_manager& dmgr)
  {
    typename MPI_detail::work_item_traits<WorkItem>::incoming_payload_type payload;
    
    this->mpi_world.recv(source, MPI_detail::MESSAGE_WORK_PRODUCT_READY, payload);
    dmgr.store(token, payload.get_data());
  }


#endif //LSSEFT_MASTER_CONTROLLER_H
