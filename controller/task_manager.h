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

#ifndef LSSEFT_TASK_MANAGER_H
#define LSSEFT_TASK_MANAGER_H


#include <memory>


#include "argument_cache.h"
#include "master_controller.h"
#include "slave_controller.h"

#include "error/error_handler.h"

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
    boost::mpi::environment mpi_env;

    //! shared BOOST MPI communicator
    boost::mpi::communicator mpi_world;


    // Caches

    //! argument cache
    argument_cache arg_cache;


    // Controller delegates
    // note these must be declared after all objects on which they depend for construction;
    // this guarantees that these objects will have been constructed by the same
    // the delegates are constructed

    //! master delegate
    master_controller master_ctrl;

    //! slave delegate
    slave_controller slave_ctrl;

  };


#endif //LSSEFT_TASK_MANAGER_H
