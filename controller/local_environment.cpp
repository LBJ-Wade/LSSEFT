//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include <cstdlib>

#include "local_environment.h"

#include "utilities/python_finder.h"


local_environment::local_environment()
  {
    // set up python path
    python_location = utilities::find_python();

    // determine if terminal supports colour output
    char* term_type_cstr = std::getenv("TERM");

    if(term_type_cstr == nullptr)
      {
        colour_output = false;
        return;
      }

    std::string term_type(term_type_cstr);

    colour_output = term_type == "xterm"
                    || term_type == "xterm-color"
                    || term_type == "xterm-256color"
                    || term_type == "screen"
                    || term_type == "linux"
                    || term_type == "cygwin";
  }
