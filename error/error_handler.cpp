//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include <iostream>

#include "error_handler.h"

#include "misc/ansi_colour_codes.h"

#include "localizations/messages.h"


error_handler::error_handler(argument_cache& ac, local_environment& le)
  : arg_cache(ac),
    local_env(le)
  {
  }


void error_handler::error(std::string msg)
  {
    bool colour = this->local_env.get_terminal_colour_support() && this->arg_cache.get_colour_output();

    if(colour) std::cout << ANSI_BOLD_RED;
    std::cout << msg << '\n';
    if(colour) std::cout << ANSI_NORMAL;
  }


void error_handler::warn(std::string msg)
  {
    bool colour = this->local_env.get_terminal_colour_support() && this->arg_cache.get_colour_output();

    if(colour) std::cout << ANSI_BOLD_MAGENTA;
    std::cout << WARNING_LABEL << " ";
    if(colour) std::cout << ANSI_NORMAL;
    std::cout << msg << '\n';
  }


void error_handler::info(std::string msg)
  {
    bool colour = this->local_env.get_terminal_colour_support() && this->arg_cache.get_colour_output();

    if(this->arg_cache.get_verbose())
      {
        std::cout << msg << '\n';
      }
  }
