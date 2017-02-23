//
// Created by David Seery on 21/11/2015.
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

#ifndef LSSEFT_WAVENUMBER_RECORD_H
#define LSSEFT_WAVENUMBER_RECORD_H


#include <memory>
#include <map>

#include "tokens.h"
#include "units/Mpc_units.h"

#include "boost/serialization/serialization.hpp"
#include "boost/serialization/map.hpp"
#include "boost/serialization/shared_ptr.hpp"


template <typename Token>
class wavenumber_record
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! build a wavenumber record
    wavenumber_record(const Mpc_units::energy& _k, const Token& tok);


    // INTERFACE

  public:

    //! deference to get k-value (not we return a copy, not a reference)
    Mpc_units::energy operator*() const { return(this->k); }

    //! get token
    const Token& get_token() const { return(this->token); }


    // INTERNAL DATA

  private:

    //! value of wavenumber
    Mpc_units::energy k;

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
wavenumber_record<Token>::wavenumber_record(const Mpc_units::energy& _k, const Token& tok)
  : k(_k),
    token(tok)
  {
  }


namespace boost
  {

    // wavenumber_record has no default constructor, and therefore we have to specialize
    // load/store methods for Boost::serialization

    namespace serialization
      {

        template <typename Archive, typename Token>
        inline void save_construct_data(Archive& ar, const wavenumber_record<Token>* t, const unsigned int file_version)
          {
          }


        template <typename Archive, typename Token>
        inline void load_construct_data(Archive& ar, wavenumber_record<Token>* t, const unsigned int file_version)
          {
            // construct empty object; values will be populated by standard deserialization
            Mpc_units::energy k(0.0);
            Token tk(0);
            ::new(t) wavenumber_record<Token>(k, tk);
          }


        // for use within a std::map we also need a specialization for std::pair< Mpc_units::energy, wavenumber_record >

        template <typename Archive, typename Token>
        inline void save_construct_data(Archive& ar, const std::pair< Mpc_units::energy, wavenumber_record<Token> >* t, const unsigned int file_version)
          {
            const Mpc_units::energy& k = *(t->second);
            unsigned int id = t->second.get_token().get_id();

            ar << boost::serialization::make_nvp("first", k);
            ar << boost::serialization::make_nvp("second", id);
          }


        template <typename Archive, typename Token>
        inline void load_construct_data(Archive& ar, std::pair< Mpc_units::energy, wavenumber_record<Token> >* t, const unsigned int file_version)
          {
            Mpc_units::energy k(0);
            unsigned int id;

            ar >> boost::serialization::make_nvp("first", k);
            ar >> boost::serialization::make_nvp("second", id);

            ::new(t) std::pair< Mpc_units::energy, wavenumber_record<Token> >(k, wavenumber_record<Token>(k, id));
          }

      }   // namespace serialization

  }   // namespace boost


#endif //LSSEFT_WAVENUMBER_RECORD_H
