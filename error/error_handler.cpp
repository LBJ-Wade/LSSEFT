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


#include <iostream>

#include "error_handler.h"

#include "misc/ansi_colour_codes.h"

#include "localizations/messages.h"

#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/posix_time/time_formatters.hpp"


error_handler::error_handler(argument_cache& ac, local_environment& le)
  : arg_cache(ac),
    local_env(le)
  {
  }


void error_handler::error(std::string msg)
  {
    bool colour = this->local_env.get_terminal_colour_support() && this->arg_cache.get_colour_output();

    boost::posix_time::ptime now = boost::posix_time::second_clock::universal_time();
    std::string nowstr = boost::posix_time::to_simple_string(now);

    if(colour) std::cout << ANSI_BOLD_RED;
    std::cout << "lsseft [" << nowstr << "]: " << msg << '\n';
    if(colour) std::cout << ANSI_NORMAL;
  }


void error_handler::warn(std::string msg)
  {
    bool colour = this->local_env.get_terminal_colour_support() && this->arg_cache.get_colour_output();

    boost::posix_time::ptime now = boost::posix_time::second_clock::universal_time();
    std::string nowstr = boost::posix_time::to_simple_string(now);

    if(colour) std::cout << ANSI_BOLD_MAGENTA;
    std::cout << WARNING_LABEL << " ";
    if(colour) std::cout << ANSI_NORMAL;
    std::cout << "lsseft [" << nowstr << "]: " << msg << '\n';
  }


void error_handler::info(std::string msg)
  {
    bool colour = this->local_env.get_terminal_colour_support() && this->arg_cache.get_colour_output();

    boost::posix_time::ptime now = boost::posix_time::second_clock::universal_time();
    std::string nowstr = boost::posix_time::to_simple_string(now);

    if(this->arg_cache.get_verbose())
      {
        std::cout << "lsseft [" << nowstr << "]: " << msg << '\n';
      }
  }
