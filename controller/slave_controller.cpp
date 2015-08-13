//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "slave_controller.h"


#include "MPI_detail/mpi_operations.h"

#include "boost/mpi.hpp"


slave_controller::slave_controller(boost::mpi::environment& me, boost::mpi::communicator& mw, argument_cache& ac)
  : mpi_env(me),
    mpi_world(mw),
    arg_cache(ac),
    local_env(),
    err_handler(arg_cache, local_env)
  {
  }


void slave_controller::execute()
  {
    // wait until we receive an EXIT instruction from master
    bool finished = false;

    while(!finished)
      {
        // wait for a message from the master node
        boost::mpi::status stat = this->mpi_world.probe(MPI_detail::RANK_MASTER);

        switch(stat.tag())
          {
            case MPI_detail::MESSAGE_TERMINATE:
              this->mpi_world.recv(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_TERMINATE);
              finished = true;
              break;

            default:
              assert(false);
          }
      }
  }
