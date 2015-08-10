//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "slave_controller.h"


slave_controller::slave_controller(std::shared_ptr<boost::mpi::environment>& me,
                                   std::shared_ptr<boost::mpi::communicator>& mw, std::shared_ptr<argument_cache>& ac)
  : mpi_env(me),
    mpi_world(mw),
    arg_cache(ac)
  {
    local_env = std::make_shared<local_environment>();
  }
