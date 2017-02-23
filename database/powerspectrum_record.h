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

#ifndef LSSEFT_POWERSPECTRUM_RECORD_H
#define LSSEFT_POWERSPECTRUM_RECORD_H


#include <memory>
#include <vector>
#include <map>

#include "units/Mpc_units.h"

#include "boost/timer/timer.hpp"
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/shared_ptr.hpp"


template <typename Dimension>
class Pk_record
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    Pk_record(const Mpc_units::energy& _k, const Dimension& _Pk);

    //! destructor is default
    ~Pk_record() = default;


    // INTERFACE

  public:

    //! dereference to get Pk-value
    const Dimension& operator*() const { return(this->Pk); }
    
    //! get Pk-value
    const Dimension& get_Pk() const { return(this->Pk); }

    //! get wavenumber
    const Mpc_units::energy& get_wavenumber() const { return(this->k); }


    // INTERNAL DATA

  private:

    //! k-value in units of eV
    Mpc_units::energy k;

    //! P(k) for this k-value in units of (Mpc/h)^3
    Dimension Pk;


    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & k;
        ar & Pk;
      }

  };


template <typename Dimension>
Pk_record<Dimension>::Pk_record(const Mpc_units::energy& _k, const Dimension& _Pk)
  : k(_k),
    Pk(_Pk)
  {
  }


namespace boost
  {

    // Pk_record has no default constructor, and therefore we have to specialize
    // load/store methods for Boost::serialization

    namespace serialization
      {

        template <typename Archive, typename Dimension>
        inline void save_construct_data(Archive& ar, const Pk_record<Dimension>* t, const unsigned int file_version)
          {
          }


        template <typename Archive, typename Dimension>
        inline void load_construct_data(Archive& ar, Pk_record<Dimension>* t, const unsigned int file_version)
          {
            // construct an empty Pk_record<Dimension>() object; its values will be populated by
            // standard deserialization
            Mpc_units::energy k(0.0);
            Dimension value(0.0);
            ::new(t) Pk_record<Dimension>(k, value);
          }


        // for use within a std::map we also need a specialization for std::pair< Mpc_units::energy, Pk_record >

        template <typename Archive, typename Dimension>
        inline void save_construct_data(Archive& ar, const std::pair< const Mpc_units::energy, Pk_record<Dimension> >* t, unsigned int file_version)
          {
            const Mpc_units::energy& k = t->second.get_wavenumber();
            const Dimension& value = t->second.get_Pk();

            ar << boost::serialization::make_nvp("first", k);
            ar << boost::serialization::make_nvp("second", value);
          }


        template <typename Archive, typename Dimension>
        inline void load_construct_data(Archive& ar, std::pair< const Mpc_units::energy, Pk_record<Dimension> >* t, unsigned int file_version)
          {
            Mpc_units::energy k(0.0);
            Dimension value(0.0);

            ar >> boost::serialization::make_nvp("first", k);
            ar >> boost::serialization::make_nvp("second", value);

            ::new(t) std::pair< Mpc_units::energy, Pk_record<Dimension> >(k, Pk_record<Dimension>(k, value));
          }

      }   // namespace serialization

  }   // namespace boost


#endif //LSSEFT_POWERSPECTRUM_RECORD_H
