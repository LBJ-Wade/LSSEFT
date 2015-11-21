//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_WAVENUMBER_RECORD_H
#define LSSEFT_WAVENUMBER_RECORD_H


#include <memory>
#include <map>

#include "tokens.h"
#include "units/eV_units.h"

#include "boost/serialization/serialization.hpp"
#include "boost/serialization/map.hpp"
#include "boost/serialization/shared_ptr.hpp"


template <typename Token>
class wavenumber_record
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! build a wavenumber record
    wavenumber_record(const eV_units::energy& _k, const Token& tok);


    // INTERFACE

  public:

    //! deference to get k-value (not we return a copy, not a reference)
    eV_units::energy operator*() const { return(this->k); }

    //! get token
    const Token& get_token() const { return(this->token); }


    // INTERNAL DATA

  private:

    //! value of wavenumber
    eV_units::energy k;

    //! database token
    Token token;


    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & k;
        ar & token;
      }

  };


template <typename Token>
wavenumber_record<Token>::wavenumber_record(const eV_units::energy& _k, const Token& tok)
  : k(_k),
    token(tok)
  {
  }


#endif //LSSEFT_WAVENUMBER_RECORD_H
