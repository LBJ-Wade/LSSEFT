//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_REDSHIFT_DATABASE_H
#define LSSEFT_REDSHIFT_DATABASE_H


#include <memory>
#include <map>

#include "tokens.h"


class redshift_record
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! build a redshift record
    redshift_record(double _z, std::shared_ptr<redshift_token> tok);


    // INTERNAL DATA

  private:

    //! value of redshift
    double z;

    //! database token
    std::shared_ptr<redshift_token> token;

  };


class redshift_database
  {

  private:

    //! alias for data structure
    typedef std::map< unsigned int, redshift_record > database_type;


    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    redshift_database() = default;

    //! destructor
    ~redshift_database() = default;


    // INTERFACE -- ADD AND LOOKUP RECORDS

  public:

    //! add record to the database

    //! The record shouldn't already exist. No checks are made to test for duplicates
    void add_record(double z, std::shared_ptr<redshift_token> tok);


    // INTERNAL DATA

  private:

    //! database
    database_type database;

  };


#endif //LSSEFT_REDSHIFT_DATABASE_H
