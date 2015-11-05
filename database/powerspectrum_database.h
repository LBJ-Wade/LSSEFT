//
// Created by David Seery on 04/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_POWERSPECTRUM_DATABASE_H
#define LSSEFT_POWERSPECTRUM_DATABASE_H

#include <memory>
#include <vector>
#include <map>

#include "generic_record_iterator.h"
#include "generic_value_iterator.h"

#include "units/eV_units.h"

#include "boost/timer/timer.hpp"
#include "boost/serialization/serialization.hpp"


class Pk_record
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    Pk_record(const eV_units::energy& in_k, double in_Pk)
      : k(in_k),
        Pk(in_Pk)
      {
      }

    //! destructor is default
    ~Pk_record() = default;


    // INTERFACE

  public:

    //! dereference to get Pk-value (note we return a copy, not a reference)
    double operator*() const { return(this->Pk); }

    //! get wavenumber
    const eV_units::energy& get_wavenumber() const { return(this->k); }

    //! get Pk-value
    double get_Pk() const { return(this->Pk); }


    // INTERNAL DATA

  private:

    eV_units::energy k;
    double Pk;

  };


class powerspectrum_database
  {

  private:

    // use a map type to store the ordered power spectrum data
    // so a power spectrum is a map k -> P(k)
    // where k is measured in eV
    typedef std::map< eV_units::energy, Pk_record > database_type;


    // RECORD-VALUED ITERATORS

  public:

    // specialize generic_record_iterator<> to obtain const and non-const iterators into the database
    typedef configuration_database::generic_record_iterator< database_type::iterator, database_type::const_iterator, Pk_record, false > record_iterator;
    typedef configuration_database::generic_record_iterator< database_type::iterator, database_type::const_iterator, Pk_record, true >  const_record_iterator;

    typedef configuration_database::generic_record_iterator< database_type::reverse_iterator, database_type::const_reverse_iterator, Pk_record, false > reverse_record_iterator;
    typedef configuration_database::generic_record_iterator< database_type::reverse_iterator, database_type::const_reverse_iterator, Pk_record, true >  const_reverse_record_iterator;


    // CONFIGURATION-VALUED ITERATORS

  public:

    typedef configuration_database::generic_value_iterator< database_type::iterator, database_type::const_iterator, double, false > value_iterator;
    typedef configuration_database::generic_value_iterator< database_type::iterator, database_type::const_iterator, double, true >  const_value_iterator;

    typedef configuration_database::generic_value_iterator< database_type::reverse_iterator, database_type::const_reverse_iterator, double, false > reverse_value_iterator;
    typedef configuration_database::generic_value_iterator< database_type::reverse_iterator, database_type::const_reverse_iterator, double, true >  const_reverse_value_iterator;


    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor is default
    powerspectrum_database() = default;

    //! destructor is default
    ~powerspectrum_database() = default;


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

    //! the record shouldn't already exist, but no check is made to enforce this
    void add_record(const eV_units::energy& k, double Pk);


    // UTILITY FUNCTIONS

  public:

    //! get number of elements in the database
    size_t size() const { return(this->database.size()); }


    // INTERNAL DATA

  private:

    //! database cache
    database_type database;

  };


#endif //LSSEFT_POWERSPECTRUM_DATABASE_H
