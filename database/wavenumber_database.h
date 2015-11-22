//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_WAVENUMBER_DATABASE_H
#define LSSEFT_WAVENUMBER_DATABASE_H


#include <memory>
#include <map>

#include "tokens.h"
#include "generic_record_iterator.h"
#include "generic_value_iterator.h"
#include "wavenumber_record.h"

#include "units/Mpc_units.h"

#include "boost/serialization/serialization.hpp"
#include "boost/serialization/map.hpp"
#include "boost/serialization/split_member.hpp"


template <typename Token>
class wavenumber_database
  {

  private:

    //! alias for data structure;
    //! records are stored in ascending wavenumber order
    typedef std::map< Mpc_units::energy, wavenumber_record<Token> > database_type;

    //! alias for lookup-by-key index
    typedef std::map< unsigned int, typename database_type::iterator > key_index_type;


    // RECORD-VALUED ITERATORS

  public:

    // specialize generic_record_iterator<> to obtain const and non-const iterators into the database
    typedef configuration_database::generic_record_iterator< typename database_type::iterator, typename database_type::const_iterator, wavenumber_record<Token>, false > record_iterator;
    typedef configuration_database::generic_record_iterator< typename database_type::iterator, typename database_type::const_iterator, wavenumber_record<Token>, true >  const_record_iterator;

    typedef configuration_database::generic_record_iterator< typename database_type::reverse_iterator, typename database_type::const_reverse_iterator, wavenumber_record<Token>, false > reverse_record_iterator;
    typedef configuration_database::generic_record_iterator< typename database_type::reverse_iterator, typename database_type::const_reverse_iterator, wavenumber_record<Token>, true >  const_reverse_record_iterator;


    // CONFIGURATION-VALUED ITERATORS

  public:

    typedef configuration_database::generic_value_iterator< typename database_type::iterator, typename database_type::const_iterator, Mpc_units::energy, false > value_iterator;
    typedef configuration_database::generic_value_iterator< typename database_type::iterator, typename database_type::const_iterator, Mpc_units::energy, true >  const_value_iterator;

    typedef configuration_database::generic_value_iterator< typename database_type::reverse_iterator, typename database_type::const_reverse_iterator, Mpc_units::energy, false > reverse_value_iterator;
    typedef configuration_database::generic_value_iterator< typename database_type::reverse_iterator, typename database_type::const_reverse_iterator, Mpc_units::energy, true >  const_reverse_value_iterator;


    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    wavenumber_database() = default;

    //! destructor
    ~wavenumber_database() = default;


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
    void add_record(const Mpc_units::energy& k, const Token& tok);

    //! lookup record by token
    record_iterator lookup(Token tok);

    //! lookup record by token -- const version
    const_record_iterator lookup(Token tok) const;

  protected:

    //! rebuild key index
    void rebuild_key_index();


    // UTILITY FUNCTIONS

  public:

    //! get number of elements in the database
    size_t size() const { return(this->database.size()); }


    // INTERNAL DATA

  private:

    //! database
    database_type database;

    //! lookup-by-key index
    key_index_type key_index;


    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void save(Archive& ar, unsigned int version) const
      {
        ar & database;
      }

    template <typename Archive>
    void load(Archive& ar, unsigned int version)
      {
        ar & database;
        this->rebuild_key_index();
      }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

  };


template <typename Token>
void wavenumber_database<Token>::add_record(const Mpc_units::energy& k, const Token& tok)
  {
    std::pair< typename database_type::iterator, bool> emplaced_value = this->database.emplace(k, wavenumber_record<Token>(k, tok));
    assert(emplaced_value.second);

    this->key_index.emplace(tok.get_id(), emplaced_value.first);
  }


template <typename Token>
typename wavenumber_database<Token>::record_iterator wavenumber_database<Token>::lookup(Token tok)
  {
    typename key_index_type::iterator t = this->key_index.find(tok.get_id());
    typename database_type::iterator dt = t->second;

    return record_iterator(dt);
  }


template <typename Token>
typename wavenumber_database<Token>::const_record_iterator wavenumber_database<Token>::lookup(Token tok) const
  {
    typename key_index_type::const_iterator t = this->key_index.find(tok.get_id());
    typename database_type::const_iterator dt = t->second;

    return const_record_iterator(dt);
  }


template <typename Token>
void wavenumber_database<Token>::rebuild_key_index()
  {
    this->key_index.clear();

    for(typename database_type::iterator t = this->database.begin(); t != this->database.end(); ++t)
      {
        this->key_index.emplace(t->second.get_token().get_id(), t);
      }
  }


#endif //LSSEFT_WAVENUMBER_DATABASE_H
