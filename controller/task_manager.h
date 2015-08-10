//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_TASK_MANAGER_H
#define LSSEFT_TASK_MANAGER_H


#include <memory>


#include "argument_cache.h"

#include "master_controller.h"
#include "slave_controller.h"


#include "boost/mpi.hpp"


class task_manager
  {

  public:

    //! construct a task manager object using command line arguments
    task_manager(int argc, char* argv[]);

    //! destructor is default (MPI environment management is done via boost::mpi::environment, whose destructor finalizes the MPI session)
    ~task_manager() = default;


    // INTERFACE

  public:

    //! process work
    void work();


    // INTERNAL DATA

  private:

    // MPI environment

    //! shared BOOST MPI environment
    std::shared_ptr<boost::mpi::environment> mpi_env;

    //! shared BOOST MPI communicator
    std::shared_ptr<boost::mpi::communicator> mpi_world;


    // Caches

    //! argument cache
    std::shared_ptr<argument_cache> arg_cache;


    // Controller delegates
    // note these must be declared after all shared_ptr<> objects on which they depend for construction;
    // this guarantees that these objects will have been constructed by the same
    // the delegates are constructed

    //! master delegate
    master_controller master_ctrl;

    //! slave delegate
    slave_controller slave_ctrl;

  };


#endif //LSSEFT_TASK_MANAGER_H
