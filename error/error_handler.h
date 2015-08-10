//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
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
    error_handler(std::shared_ptr<argument_cache>& ac,
                  std::shared_ptr<local_environment>& le);

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
    std::shared_ptr<argument_cache> arg_cache;

    //! local environment object
    std::shared_ptr<local_environment> local_env;

  };


#endif //LSSEFT_ERROR_HANDLER_H
