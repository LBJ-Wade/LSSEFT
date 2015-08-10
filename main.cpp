//
// Created by David Seery on 31/07/15.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include <iostream>

#include "controller/task_manager.h"

int main(int argc, char* argv[])
  {
    // instantiate task manager
    task_manager mgr(argc, argv);

    mgr.work();

    return(EXIT_SUCCESS);
  }
