//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "task_manager.h"

#include "MPI_detail/mpi_operations.h"


//! construct a task manager object
task_manager::task_manager(int argc, char* argv[])
  : mpi_env(argc, argv),
    mpi_world(),
    arg_cache(),
    master_ctrl(mpi_env, mpi_world, arg_cache),
    slave_ctrl(mpi_env, mpi_world, arg_cache)
  {
    if(mpi_world.rank() == MPI_detail::RANK_MASTER)
      {
        master_ctrl.process_arguments(argc, argv);
      }
  }


void task_manager::work()
  {
    if(this->mpi_world.rank() == MPI_detail::RANK_MASTER)
      {
        this->master_ctrl.execute();
      }
    else
      {
        this->slave_ctrl.execute();
      }
  }
