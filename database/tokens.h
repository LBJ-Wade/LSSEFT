//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_TOKENS_H
#define LSSEFT_TOKENS_H


//! base class: a generic token
class generic_token
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    generic_token(unsigned int i);

    //! destructor is default
    virtual ~generic_token() = default;


    // INTERFACE

  public:

    virtual unsigned int get_id() const { return(this->id); }


    // INTERNAL DATA

  protected:

    unsigned int id;

  };


//! token representing an FRW model
class FRW_model_token: public generic_token
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    FRW_model_token(unsigned int i);

    //! destructor is default
    virtual ~FRW_model_token() = default;

  };


#endif //LSSEFT_TOKENS_H
