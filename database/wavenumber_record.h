//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
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
            ar << *(*t);                    // store wavenumber value
            ar << t->get_token().get_id();  // store token identifier
          }


        template <typename Archive, typename Token>
        inline void load_construct_data(Archive& ar, wavenumber_record<Token>* t, const unsigned int file_version)
          {
            double k;
            unsigned int id;

            ar >> k;    // unpack wavenumber value
            ar >> id;   // unpack token identifier

            Mpc_units::energy k_in_eV(k);

            // invoke in-place constructor

            ::new(t) wavenumber_record<Token>(k_in_eV, Token(id));
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
