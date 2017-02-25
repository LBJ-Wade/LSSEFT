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

#ifndef LSSEFT_REDSHIFT_RECORD_H
#define LSSEFT_REDSHIFT_RECORD_H


#include <memory>
#include <map>

#include "tokens.h"

#include "boost/serialization/serialization.hpp"
#include "boost/serialization/map.hpp"
#include "boost/serialization/shared_ptr.hpp"


class z_record
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! build a redshift record
    z_record(double _z, const z_token& tok);


    // INTERFACE

  public:

    //! deference to get z-value (note we return a copy, not a reference)
    double operator*() const { return(this->z); }

    //! get token
    const z_token& get_token() const { return(this->token); }


    // INTERNAL DATA

  private:

    //! value of redshift
    double z;

    //! database token
    z_token token;


    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & z;
        ar & token;
      }

  };


namespace boost
  {

    // z_record has no default constructor, and therefore we have to specialize
    // load/store methods for Boost:serialization

    namespace serialization
      {

        template <typename Archive>
        inline void save_construct_data(Archive& ar, const z_record* t, const unsigned int file_version)
          {
          }


        template <typename Archive>
        inline void load_construct_data(Archive& ar, z_record* t, const unsigned int file_version)
          {
            // construct empty object; values will be populated later by standard deserialization
            z_token tk(0);
            ::new(t) z_record(0.0, tk);
          }


        // for use within a std::map (eg in z_database), we also need a specialization for std::pair< unsigned int, z_record >

        template <typename Archive>
        inline void save_construct_data(Archive& ar, const std::pair< const double, z_record>* t, const unsigned int file_version)
          {
            double z = *(t->second);
            unsigned int id = t->second.get_token().get_id();

            ar << boost::serialization::make_nvp("first", z);
            ar << boost::serialization::make_nvp("second", id);
          }


        template <typename Archive>
        inline void load_construct_data(Archive& ar, std::pair< const double, z_record>* t, const unsigned int file_version)
          {
            double z;
            unsigned int id;

            ar >> boost::serialization::make_nvp("first", z);
            ar >> boost::serialization::make_nvp("second", id);

            ::new(t) std::pair< const double, z_record>(z, z_record(z, id));
          }

      }   // namespace serialization

  }   // namespace boost


#endif //LSSEFT_REDSHIFT_RECORD_H
