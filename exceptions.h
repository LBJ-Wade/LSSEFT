//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_EXCEPTIONS_H
#define LSSEFT_EXCEPTIONS_H


#include <stdexcept>
#include <string>


enum class exception_type
  {
    database_error,
    sqlite3_error,
    transaction_error
  };


class runtime_exception: public std::runtime_error
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    runtime_exception(exception_type t, const std::string msg)
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
