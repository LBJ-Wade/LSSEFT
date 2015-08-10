//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_MASTER_CONTROLLER_H
#define LSSEFT_MASTER_CONTROLLER_H


#include <memory>

#include "argument_cache.h"
#include "local_environment.h"

#include "error/error_handler.h"

#include "boost/mpi.hpp"


class master_controller
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! construct a master controller delegate
    master_controller(std::shared_ptr<boost::mpi::environment>& me,
                      std::shared_ptr<boost::mpi::communicator>& mw,
                      std::shared_ptr<argument_cache>& ac);

    //! destructor is default
    ~master_controller() = default;


    // INTERFACE

  public:

    //! process command-line arguments
    void process_arguments(int argc, char* argv[]);

    //! execute
    void execute();


    // INTERNAL DATA

  private:

    // MPI environment

    //! MPI environment object, inherited from parent task manager
    std::shared_ptr<boost::mpi::environment> mpi_env;

    //! MPI communicator, inherited from parent task manager
    std::shared_ptr<boost::mpi::communicator> mpi_world;


    // Caches

    //! argument cache, inherited from parent task manager
    std::shared_ptr<argument_cache> arg_cache;

    //! local environment properties (constructed locally)
    std::shared_ptr<local_environment> local_env;


    // Functional blocks

    //! error handler
    std::shared_ptr<error_handler> err_handler;

  };


#endif //LSSEFT_MASTER_CONTROLLER_H
