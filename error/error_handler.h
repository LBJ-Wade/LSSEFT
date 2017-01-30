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

#ifndef LSSEFT_ERROR_HANDLER_H
#define LSSEFT_ERROR_HANDLER_H


#include <memory>

#include "controller/argument_cache.h"
#include "controller/local_environment.h"


class error_handler
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    error_handler(argument_cache& ac, local_environment& le);

    //! destructor is default
    ~error_handler() = default;


    // INTERFACE

  public:

    //! report an error
    void error(std::string msg);

    //! report a warning
    void warn(std::string msg);

    //! report an information message
    void info(std::string msg);


    // INTERNAL DATA

  private:

    //! argument cache inherited from parent controller
    argument_cache& arg_cache;

    //! local environment object
    local_environment& local_env;

  };


#endif //LSSEFT_ERROR_HANDLER_H
