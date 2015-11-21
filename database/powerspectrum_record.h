//
// Created by David Seery on 21/11/2015.
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

  };


#endif //LSSEFT_POWERSPECTRUM_RECORD_H
