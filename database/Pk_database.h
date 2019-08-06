//
// Created by David Seery on 04/11/2015.
// --@@ // Copyright (c) 2017 University of Sussex. All rights reserved.
//
// This file is part of the Sussex Effective Field Theory for
// Large-Scale Structure platform (LSSEFT).
//
// LSSEFT is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// LSSEFT is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with LSSEFT.  If not, see <http://www.gnu.org/licenses/>.
//
// @license: GPL-2
// @contributor: David Seery <D.Seery@sussex.ac.uk>
// --@@
//

#ifndef LSSEFT_POWERSPECTRUM_DATABASE_H
#define LSSEFT_POWERSPECTRUM_DATABASE_H

#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <map>

#include "generic_record_iterator.h"
#include "generic_value_iterator.h"
#include "powerspectrum_record.h"

#include "units/Mpc_units.h"

#include "exceptions.h"
#include "localizations/messages.h"

#include "boost/timer/timer.hpp"
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/map.hpp"
#include "boost/filesystem.hpp"


namespace Pk_database_impl
  {
    
    
    template<typename Dimension> struct DimensionTraits;
    
    template <>
    struct DimensionTraits<Mpc_units::inverse_energy3>
      {
        constexpr Mpc_units::inverse_energy3 unit() const { return Mpc_units::Mpc3; }
      };
    
    template <>
    struct DimensionTraits<Mpc_units::inverse_energy>
      {
        constexpr Mpc_units::inverse_energy unit() const { return Mpc_units::Mpc; }
      };
    
    
  }   // namespace Pk_database_impl


template <typename Dimension>
class Pk_database
  {

  private:

    // use a map type to store the ordered power spectrum data
    // so a power spectrum is a map k -> P(k)
    // where k is measured in eV
    typedef std::map< Mpc_units::energy, Pk_record<Dimension> > database_type;

    typedef Dimension Pk_units;


    // RECORD-VALUED ITERATORS

  public:

    // specialize generic_record_iterator<> to obtain const and non-const iterators into the database
    typedef configuration_database::generic_record_iterator< typename database_type::iterator, typename database_type::const_iterator, Pk_record<Dimension>, false > record_iterator;
    typedef configuration_database::generic_record_iterator< typename database_type::iterator, typename database_type::const_iterator, Pk_record<Dimension>, true >  const_record_iterator;

    typedef configuration_database::generic_record_iterator< typename database_type::reverse_iterator, typename database_type::const_reverse_iterator, Pk_record<Dimension>, false > reverse_record_iterator;
    typedef configuration_database::generic_record_iterator< typename database_type::reverse_iterator, typename database_type::const_reverse_iterator, Pk_record<Dimension>, true >  const_reverse_record_iterator;


    // CONFIGURATION-VALUED ITERATORS

  public:

    typedef configuration_database::generic_value_iterator< typename database_type::iterator, typename database_type::const_iterator, Pk_units, false > value_iterator;
    typedef configuration_database::generic_value_iterator< typename database_type::iterator, typename database_type::const_iterator, Pk_units, true >  const_value_iterator;

    typedef configuration_database::generic_value_iterator< typename database_type::reverse_iterator, typename database_type::const_reverse_iterator, Pk_units, false > reverse_value_iterator;
    typedef configuration_database::generic_value_iterator< typename database_type::reverse_iterator, typename database_type::const_reverse_iterator, Pk_units, true >  const_reverse_value_iterator;


    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! empty constructor
    Pk_database();
    
    //! construct a power spectrum database from a file in CAMB format
    Pk_database(const boost::filesystem::path& p);

    //! destructor is default
    ~Pk_database() = default;


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
    void add_record(const Mpc_units::energy& k, const Dimension& Pk);


    // UTILITY FUNCTIONS

  public:

    //! get number of elements in the database
    size_t size() const { return(this->database.size()); }

    //! get largest stored k-value
    const Mpc_units::energy& get_k_min() const { return(this->k_min); }

    //! get smallest stored k-value
    const Mpc_units::energy& get_k_max() const { return(this->k_max); }


    // INTERNAL DATA

  private:

    //! database cache
    database_type database;

    //! smallest k-mode in the database
    Mpc_units::energy k_min;

    //! largest k-mode in the database
    Mpc_units::energy k_max;


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


template <typename Dimension>
Pk_database<Dimension>::Pk_database()
  : k_min(std::numeric_limits<double>::max()),
    k_max(std::numeric_limits<double>::min())
  {
  }


template <typename Dimension>
Pk_database<Dimension>::Pk_database(const boost::filesystem::path& p)
  : Pk_database<Dimension>()   // forward to empty constructor
  {
    // set up an input stream to ingest the file, and open it
    std::ifstream in;
    in.open(p.string());

    // check whether the input stream is in good condition
    if(!in.good())
      {
        std::ostringstream msg;
        msg << ERROR_POWERSPECTRUM_FILE_NOT_READABLE_A << " " << p << " " << ERROR_POWERSPECTRUM_FILE_NOT_READABLE_B;
        throw runtime_exception(exception_type::runtime_error, msg.str());
      }

    unsigned line_number = 0;

    for(std::string line; std::getline(in, line);)
      {
        ++line_number;

        if(!line.empty())
          {
            std::stringstream line_stream(line);

            if(line.front() != '#')   // hash # is CAMB-format comment character
              {
                double _k, _Pk;
                line_stream >> _k >> _Pk;

                // check whether extraction was successful
                if(!line_stream.fail())
                  {
                    Mpc_units::energy k = _k / Mpc_units::Mpc;
                    Dimension Pk = _Pk * Pk_database_impl::DimensionTraits<Dimension>().unit();
                    this->add_record(k, Pk);
                  }
                else
                  {
                    std::cerr << "Ignored from input power spectrum '" << p.string() << "'\n";
                    std::cerr << "  " << line_number << ": " << line << "\n";
                  }
              }
          }
      }

    in.close();
  }


template <typename Dimension>
void Pk_database<Dimension>::add_record(const Mpc_units::energy& k, const Dimension& Pk)
  {
    std::pair<typename database_type::iterator, bool> emplaced_value = this->database.emplace(k, Pk_record<Dimension>(k, Pk));
    assert(emplaced_value.second);
    
    if(k > this->k_max) this->k_max = k;
    if(k < this->k_min) this->k_min = k;
  }


#endif //LSSEFT_POWERSPECTRUM_DATABASE_H
