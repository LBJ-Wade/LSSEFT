//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SLAVE_CONTROLLER_H
#define LSSEFT_SLAVE_CONTROLLER_H


#include <memory>

#include "argument_cache.h"
#include "local_environment.h"

#include "MPI_detail/mpi_operations.h"

#include "error/error_handler.h"

#include "boost/mpi.hpp"


class slave_controller
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! construct a master controller delegate
    slave_controller(boost::mpi::environment& me, boost::mpi::communicator& mw, argument_cache& ac);

    //! destructor is default
    ~slave_controller() = default;


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


    // INTERFACE

  public:

    //! execute
    void execute();


    // TASK HANDLING

  protected:

    //! process task corresponding to given work item type
    template <typename WorkItem>
    void process_task();


    // TRANSFER FUNCTION TASKS

  protected:

    //! integrate a given transfer function
    void process_item(MPI_detail::new_transfer_integration& payload);


    // LOOP MOMENTUM TASKS

  protected:

    //! integrate a given loop
    void process_item(MPI_detail::new_loop_momentum_integration& payload);
    
    //! combine loop integral and growth-factor data to produce a 1-loop power spectrum
    void process_item(MPI_detail::new_one_loop_Pk& payload);
    
    //! combine 1-loop power spectrum data to produce multipole power spectra
    void process_item(MPI_detail::new_multipole_Pk& payload);


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


#endif //LSSEFT_SLAVE_CONTROLLER_H
