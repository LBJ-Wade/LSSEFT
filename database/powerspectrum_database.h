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
#include "powerspectrum_record.h"

#include "units/eV_units.h"

#include "boost/timer/timer.hpp"
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/map.hpp"


class powerspectrum_database
  {

  private:

    // use a map type to store the ordered power spectrum data
    // so a power spectrum is a map k -> P(k)
    // where k is measured in eV
    typedef std::map< eV_units::energy, Pk_record > database_type;

    typedef eV_units::inverse_energy3 Pk_units;


    // RECORD-VALUED ITERATORS

  public:

    // specialize generic_record_iterator<> to obtain const and non-const iterators into the database
    typedef configuration_database::generic_record_iterator< database_type::iterator, database_type::const_iterator, Pk_record, false > record_iterator;
    typedef configuration_database::generic_record_iterator< database_type::iterator, database_type::const_iterator, Pk_record, true >  const_record_iterator;

    typedef configuration_database::generic_record_iterator< database_type::reverse_iterator, database_type::const_reverse_iterator, Pk_record, false > reverse_record_iterator;
    typedef configuration_database::generic_record_iterator< database_type::reverse_iterator, database_type::const_reverse_iterator, Pk_record, true >  const_reverse_record_iterator;


    // CONFIGURATION-VALUED ITERATORS

  public:

    typedef configuration_database::generic_value_iterator< database_type::iterator, database_type::const_iterator, Pk_units, false > value_iterator;
    typedef configuration_database::generic_value_iterator< database_type::iterator, database_type::const_iterator, Pk_units, true >  const_value_iterator;

    typedef configuration_database::generic_value_iterator< database_type::reverse_iterator, database_type::const_reverse_iterator, Pk_units, false > reverse_value_iterator;
    typedef configuration_database::generic_value_iterator< database_type::reverse_iterator, database_type::const_reverse_iterator, Pk_units, true >  const_reverse_value_iterator;


    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    powerspectrum_database();

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
    void add_record(const eV_units::energy& k, const eV_units::inverse_energy3& Pk);


    // UTILITY FUNCTIONS

  public:

    //! get number of elements in the database
    size_t size() const { return(this->database.size()); }

    //! get largest stored k-value
    const eV_units::energy& get_k_min() const { return(this->k_min); }

    //! get smallest stored k-value
    const eV_units::energy& get_k_max() const { return(this->k_max); }


    // INTERNAL DATA

  private:

    //! database cache
    database_type database;

    //! smallest k-mode in the database
    eV_units::energy k_min;

    //! largest k-mode in the database
    eV_units::energy k_max;


    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & database;
        ar & k_min;
        ar & k_max;
      }

  };


#endif //LSSEFT_POWERSPECTRUM_DATABASE_H
