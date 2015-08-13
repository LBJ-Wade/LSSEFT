//
// Created by David Seery on 10/08/2015.
//

#ifndef LSSEFT_SLAVE_CONTROLLER_H
#define LSSEFT_SLAVE_CONTROLLER_H


#include <memory>

#include "argument_cache.h"
#include "local_environment.h"

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
