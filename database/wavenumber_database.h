//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_WAVENUMBER_DATABASE_H
#define LSSEFT_WAVENUMBER_DATABASE_H


#include <memory>
#include <map>

#include "tokens.h"

#include "units/eV_units.h"


class wavenumber_record
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! build a wavenumber record
    wavenumber_record(const eV_units::energy& _k, std::shared_ptr<wavenumber_token> tok);


    // INTERNAL DATA

  private:

    //! value of wavenumber
    eV_units::energy k;

    //! database token
    std::shared_ptr<wavenumber_token> token;

  };


class wavenumber_database
  {

  private:

    //! alias for data structure
    typedef std::map< unsigned int, wavenumber_record > database_type;

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    wavenumber_database() = default;

    //! destructor
    ~wavenumber_database() = default;


    // INTERFACE -- ADD AND LOOKUP RECORDS

  public:

    //! add record to the database

    //! The record shouldn't already exist. No checks are made to test for duplicates
    void add_record(const eV_units::energy& k, std::shared_ptr<wavenumber_token> tok);


    // INTERNAL DATA

  private:

    //! database
    database_type database;

  };


#endif //LSSEFT_WAVENUMBER_DATABASE_H
