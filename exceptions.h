//
// Created by David Seery on 11/08/2015.
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

#ifndef LSSEFT_EXCEPTIONS_H
#define LSSEFT_EXCEPTIONS_H


#include <stdexcept>
#include <string>


enum class exception_type
  {
    database_error,
    sqlite3_error,
    transaction_error,
    runtime_error,
    spline_error,
    filter_failure,
    store_error
  };


class runtime_exception: public std::runtime_error
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor with message
    runtime_exception(exception_type t, const std::string msg="")
      : type(t),
        std::runtime_error(msg)
      {
      }

    //! destructor is default
    virtual ~runtime_exception() = default;


    // INTERFACE

  public:

    exception_type get_exception_code() { return(this->type); }


    // INTERNAL DATA

  protected:

    exception_type type;

  };


#endif //LSSEFT_EXCEPTIONS_H
