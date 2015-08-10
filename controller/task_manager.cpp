//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "task_manager.h"

#include "MPI/mpi_operations.h"


//! construct a task manager object
//! first, all shared_ptr<> objects must be declared followed by the controller delegates
task_manager::task_manager(int argc, char* argv[])
  : mpi_env(std::make_shared<boost::mpi::environment>(argc, argv)),
    mpi_world(std::make_shared<boost::mpi::communicator>()),
    arg_cache(std::make_shared<argument_cache>()),
    master_ctrl(mpi_env, mpi_world, arg_cache),
    slave_ctrl(mpi_env, mpi_world, arg_cache)
  {
    if(mpi_world->rank() == MPI::RANK_MASTER)
      {
        master_ctrl.process_arguments(argc, argv);
      }
  }


void task_manager::work()
  {
    if(this->mpi_world->rank() == MPI::RANK_MASTER)
      {
        std::cout << "Executing as Master" << '\n';
      }
    else
      {
        std::cout << "Executing as Slave number " << this->mpi_world->rank() << '\n';
      }
  }
