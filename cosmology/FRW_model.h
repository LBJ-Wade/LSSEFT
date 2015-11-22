//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_FRW_MODEL_H
#define LSSEFT_FRW_MODEL_H


#include "Planck_defaults.h"

#include "boost/serialization/serialization.hpp"


class FRW_model
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    //! cosmological parameters refer to their values today; if none are provided, then we use
    //! Planck values as sensible defaults
    FRW_model(double om= Planck2013::omega_m, double occ= Planck2013::omega_cc, double h_= Planck2013::h, Mpc_units::energy tc= Planck2013::T_CMB);

    //! destructor is default
    ~FRW_model() = default;


    // INTERROGATE COSMOLOGICAL PARAMETERS

  public:

    //! get value of omega_m
    double get_omega_m() const { return(this->omega_m); }

    //! get value of omega_cc
    double get_omega_cc() const { return(this->omega_cc); }

    //! get value of h
    double get_h() const { return(this->h); }

    //! get value of T_CMB
    Mpc_units::energy get_T_CMB() const { return(this->T_CMB); }


    // INTERNAL DATA

  private:

    //! omega_m, measured today
    double omega_m;

    //! omega_Lambda or omega_cc, measured today
    double omega_cc;

    //! reduced Hubble constant, ie. quoted in units of 100 km/s/Mpc, measured today
    double h;

    //! CMB temperature T_CMB, measured today
    Mpc_units::energy T_CMB;


    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & omega_m;
        ar & omega_cc;
        ar & h;
        ar & T_CMB;
      }

  };


#endif //LSSEFT_FRW_MODEL_H
