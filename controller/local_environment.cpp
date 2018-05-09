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


#include <cstdlib>

#include "defaults.h"
#include "local_environment.h"


local_environment::local_environment()
  {
    // add $PATH to object finder
    find.add_environment_variable(LSSEFT_SHELL_PATH_ENV);

    // set up python path
    auto py = find.find(LSSEFT_PYTHON_EXECUTABLE);

    if(!py)
      {
        python_location = LSSEFT_DEFAULT_PYTHON_PATH;
      }
    else
      {
        python_location = *py;
      }

    // determine if terminal supports colour output
    char* term_type_cstr = std::getenv(LSSEFT_TERM_ENV);

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
