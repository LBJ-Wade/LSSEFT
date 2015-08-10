//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
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
