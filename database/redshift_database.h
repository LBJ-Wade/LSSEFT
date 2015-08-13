//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_REDSHIFT_DATABASE_H
#define LSSEFT_REDSHIFT_DATABASE_H


#include <memory>
#include <map>

#include "tokens.h"
#include "generic_record_iterator.h"
#include "generic_value_iterator.h"


class redshift_record
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! build a redshift record
    redshift_record(double _z, std::shared_ptr<redshift_token> tok);


    // INTERFACE

  public:

    //! deference to get z-value (not we return a copy, not a reference)
    double operator*() const { return(this->z); }

    //! get token
    const std::shared_ptr<redshift_token>& get_token() const { return(this->token); }


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


    // RECORD-VALUED ITERATORS

  public:

    // specialize generic_record_iterator<> to obtain const and non-const iterators into the database
    typedef configuration_database::generic_record_iterator< database_type::iterator, database_type::const_iterator, redshift_record, false > record_iterator;
    typedef configuration_database::generic_record_iterator< database_type::iterator, database_type::const_iterator, redshift_record, true >  const_record_iterator;

    typedef configuration_database::generic_record_iterator< database_type::reverse_iterator, database_type::const_reverse_iterator, redshift_record, false > reverse_record_iterator;
    typedef configuration_database::generic_record_iterator< database_type::reverse_iterator, database_type::const_reverse_iterator, redshift_record, true >  const_reverse_record_iterator;


    // CONFIGURATION-VALUED ITERATORS

  public:

    typedef configuration_database::generic_value_iterator< database_type::iterator, database_type::const_iterator, double, false > value_iterator;
    typedef configuration_database::generic_value_iterator< database_type::iterator, database_type::const_iterator, double, true >  const_value_iterator;

    typedef configuration_database::generic_value_iterator< database_type::reverse_iterator, database_type::const_reverse_iterator, double, false > reverse_value_iterator;
    typedef configuration_database::generic_value_iterator< database_type::reverse_iterator, database_type::const_reverse_iterator, double, true >  const_reverse_value_iterator;


    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    redshift_database() = default;

    //! destructor
    ~redshift_database() = default;


    // MANUFACTURE RECORD-VALUED ITERATORS

  public:

    record_iterator               record_begin()         { return record_iterator(this->database.begin()); }
    record_iterator               record_end()           { return record_iterator(this->database.end()); }
    const_record_iterator         record_begin()   const { return const_record_iterator(this->database.cbegin()); }
    const_record_iterator         record_end()     const { return const_record_iterator(this->database.cend()); }

    const_record_iterator         record_cbegin()  const { return const_record_iterator(this->database.cbegin()); }
    const_record_iterator         record_cend()    const { return const_record_iterator(this->database.cend()); }

    reverse_record_iterator       record_rbegin()        { return reverse_record_iterator(this->database.rbegin()); }
    reverse_record_iterator       record_rend()          { return reverse_record_iterator(this->database.rend()); }
    const_reverse_record_iterator record_rbegin()  const { return const_reverse_record_iterator(this->database.crbegin()); }
    const_reverse_record_iterator record_rend()    const { return const_reverse_record_iterator(this->database.crend()); }

    const_reverse_record_iterator record_crbegin() const { return const_reverse_record_iterator(this->database.crbegin()); }
    const_reverse_record_iterator record_crend()   const { return const_reverse_record_iterator(this->database.crend()); }


    // MANUFACTURE CONFIGURATION-VALUED ITERATORS

  public:

    value_iterator               value_begin()         { return value_iterator(this->database.begin()); }
    value_iterator               value_end()           { return value_iterator(this->database.end()); }
    const_value_iterator         value_begin()   const { return const_value_iterator(this->database.cbegin()); }
    const_value_iterator         value_end()     const { return const_value_iterator(this->database.cend()); }

    const_value_iterator         value_cbegin()  const { return const_value_iterator(this->database.cbegin()); }
    const_value_iterator         value_cend()    const { return const_value_iterator(this->database.cend()); }

    reverse_value_iterator       value_rbegin()        { return reverse_value_iterator(this->database.rbegin()); }
    reverse_value_iterator       value_rend()          { return reverse_value_iterator(this->database.rend()); }
    const_reverse_value_iterator value_rbegin()  const { return const_reverse_value_iterator(this->database.crbegin()); }
    const_reverse_value_iterator value_rend()    const { return const_reverse_value_iterator(this->database.crend()); }

    const_reverse_value_iterator value_crbegin() const { return const_reverse_value_iterator(this->database.crbegin()); }
    const_reverse_value_iterator value_crend()   const { return const_reverse_value_iterator(this->database.crend()); }


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
