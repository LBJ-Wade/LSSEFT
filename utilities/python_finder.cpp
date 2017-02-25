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


#include <string>
#include <cstdio>
#include <cstdlib>

#include "python_finder.h"
#include "defaults.h"


#include "boost/algorithm/string.hpp"


namespace utilities
  {

    std::string find_python()
      {
        std::string path;
        FILE* f = popen("which python", "r");

        if(!f)
          {
            path = LSSEFT_DEFAULT_PYTHON_PATH;
          }
        else
          {
            char buffer[1024];
            char* line = fgets(buffer, sizeof(buffer), f);
            pclose(f);

            if(line != nullptr)
              {
                path = std::string(line);
                boost::algorithm::trim_right(path);
              }
            else
              {
                path = LSSEFT_DEFAULT_PYTHON_PATH;
              }
          }

        return(path);
      }

  }   // namespace utilities
