//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_POWERSPECTRUM_RECORD_H
#define LSSEFT_POWERSPECTRUM_RECORD_H


#include <memory>
#include <vector>
#include <map>

#include "units/eV_units.h"

#include "boost/timer/timer.hpp"
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/shared_ptr.hpp"


class Pk_record
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    Pk_record(const eV_units::energy& _k, const eV_units::inverse_energy3& _Pk);

    //! destructor is default
    ~Pk_record() = default;


    // INTERFACE

  public:

    //! dereference to get Pk-value (note we return a copy, not a reference)
    const eV_units::inverse_energy3& operator*() const { return(this->Pk); }

    //! get wavenumber
    const eV_units::energy& get_wavenumber() const { return(this->k); }

    //! get Pk-value
    const eV_units::inverse_energy3& get_Pk() const { return(this->Pk); }


    // INTERNAL DATA

  private:

    //! k-value in units of eV
    eV_units::energy k;

    //! P(k) for this k-value in units of (Mpc/h)^3
    eV_units::inverse_energy3 Pk;


    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
      }

  };


namespace boost
  {

    // Pk_record has no default constructor, and therefore we have to specialize
    // load/store methods for Boost::serialization

    namespace serialization
      {

        template <typename Archive>
        inline void save_construct_data(Archive& ar, const Pk_record* t, const unsigned int file_version)
          {
            const eV_units::energy& k = t->get_wavenumber();
            const eV_units::inverse_energy3& value = t->get_Pk();

            ar << k;
            ar << value;
          }


        template <typename Archive>
        inline void load_construct_data(Archive& ar, Pk_record* t, const unsigned int file_version)
          {
            eV_units::energy k(0);
            eV_units::inverse_energy3 value(0);

            ar >> k;
            ar >> value;
          }


        // for use within a std::map we also need a specialization for std::pair< eV_units::energy, Pk_record >

        template <typename Archive>
        inline void save_construct_data(Archive& ar, const std::pair< const eV_units::energy, Pk_record >* t, unsigned int file_version)
          {
            const eV_units::energy& k = t->second.get_wavenumber();
            const eV_units::inverse_energy3& value = t->second.get_Pk();

            ar << boost::serialization::make_nvp("first", k);
            ar << boost::serialization::make_nvp("second", value);
          }


        template <typename Archive>
        inline void load_construct_data(Archive& ar, std::pair< const eV_units::energy, Pk_record >* t, unsigned int file_version)
          {
            eV_units::energy k(0);
            eV_units::inverse_energy3 value(0);

            ar >> boost::serialization::make_nvp("first", k);
            ar >> boost::serialization::make_nvp("second", value);

            ::new(t) std::pair< eV_units::energy, Pk_record >(k, Pk_record(k, value));
          }

      }   // namespace serialization

  }   // namespace boost


#endif //LSSEFT_POWERSPECTRUM_RECORD_H
