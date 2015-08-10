//
// Created by David Seery on 31/07/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "argument_cache.h"


argument_cache::argument_cache()
  {
    verbose = false;
    colour_output = true;

    // no default database
    database.clear();
  }
