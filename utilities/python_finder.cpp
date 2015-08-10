//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
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
