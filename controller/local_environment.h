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

#ifndef LSSEFT_LOCAL_ENVIRONMENT_H
#define LSSEFT_LOCAL_ENVIRONMENT_H


#include "boost/filesystem/operations.hpp"


class local_environment
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor detects properties of environment
    local_environment();

    //! destructor is default
    ~local_environment() = default;


    // LOCATION OF EXECUTABLES

  public:

    std::string get_python_location() const { return(this->python_location.string()); }


    // TERMINAL PROPERTIES

  public:

    bool get_terminal_colour_support() const { return(this->colour_output); }


    // INTERNAL DATA

  protected:

    // LOCATION OF EXECUTABLES

    //! Python executable
    boost::filesystem::path python_location;


    // TERMINAL PROPERTIES

    //! terminal supports colour output?
    bool colour_output;

  };


#endif //LSSEFT_LOCAL_ENVIRONMENT_H
