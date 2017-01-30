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
